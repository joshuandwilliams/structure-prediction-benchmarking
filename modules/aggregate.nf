/*
 * =============================================================================
 * Aggregate module
 * =============================================================================
 * Final collation step.  Takes the pre-renamed per-model tagged outputs
 * from COMPUTE_METRICS (*_metrics.csv and *_best directories) and produces:
 *
 *   all_metrics.csv    — concatenated metrics across every model
 *   best_models/       — collected best PDB/CIF files from each model,
 *                        organised into <model_tag>_best/ subdirectories
 *
 * SLURM resource usage tracking is handled natively by Nextflow's trace
 * file (configured in nextflow.config) — we no longer need a custom
 * collect_slurm_stats.sh step.  Per-task elapsed/cpu/rss/vmem appear in
 * ${outdir}/trace.txt and pipeline_report.html after the run finishes.
 *
 * This process works purely off staged inputs (no reads from publishDir),
 * so -resume behaves correctly.
 */

/*
 * Inputs: tagged_metrics   — list of per-model '<tag>_metrics.csv' files.
 * All filenames are unique because each alias used a distinct model_tag, so
 * staging collisions are impossible. tagged_best_dirs — list of per-model
 * '<tag>_best' directories. Same uniqueness guarantee.
 */
process AGGREGATE_RESULTS {
    tag "${params.project_name}"
    label 'cpu'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tagged_metrics,    stageAs: 'metrics_in/*'
    path tagged_best_dirs,  stageAs: 'best_in/*'

    output:
    path "all_metrics.csv",                          emit: combined_metrics
    path "all_metrics_ranked_by_effector_rmsd.csv",  emit: ranked_metrics
    path "best_models/",                             emit: best_models_dir

    script:
    """
    set -euo pipefail

    # Helps diagnose staging issues without needing to SSH to the HPC —
    # the full listing lands in .command.err and gets echoed back by
    # Nextflow on failure.
    echo "=== Work dir contents ===" >&2
    ls -la >&2 || true
    echo "" >&2
    echo "=== metrics_in/ contents ===" >&2
    ls -la metrics_in/ 2>/dev/null >&2 || echo "(metrics_in/ missing)" >&2
    echo "" >&2
    echo "=== best_in/ contents ===" >&2
    ls -la best_in/ 2>/dev/null >&2 || echo "(best_in/ missing)" >&2
    echo "" >&2

    # Recursive find so we catch the files wherever Nextflow staged them
    # (top level, metrics_in/, or variants).  Use -L to follow the symlinks
    # that Nextflow uses for staging on NFS.  No -type filter so both
    # regular files and symlinks match.
    mapfile -t csv_files < <(find -L . -name '*_metrics.csv' ! -path '*/.command*' | sort -u)

    if [ \${#csv_files[@]} -eq 0 ]; then
        echo "ERROR: No *_metrics.csv files staged." >&2
        echo "Searched from: \$(pwd)" >&2
        find -L . ! -path '*/.command*' >&2 || true
        exit 1
    fi

    echo "Found \${#csv_files[@]} metrics CSV(s):" >&2
    printf '  %s\\n' "\${csv_files[@]}" >&2
    echo "" >&2

    first=1
    for csv in "\${csv_files[@]}"; do
        if [ "\${first}" -eq 1 ]; then
            cat "\${csv}" > all_metrics.csv
            first=0
        else
            tail -n +2 "\${csv}" >> all_metrics.csv
        fi
    done

    echo "Combined \${#csv_files[@]} metrics CSVs into all_metrics.csv"
    echo ""

    # Each staged *_best directory has a unique name.  Find them wherever
    # Nextflow placed them (top level or best_in/) and copy into
    # best_models/<model_tag>_best/.  -L follows the staging symlinks.
    mkdir -p best_models

    echo "=== Collecting best models ==="
    mapfile -t best_dirs < <(find -L . -type d -name '*_best' ! -path './best_models/*' | sort -u)

    if [ \${#best_dirs[@]} -eq 0 ]; then
        echo "  WARNING: No *_best directories found — best_models/ will be empty." >&2
    fi

    for src in "\${best_dirs[@]}"; do
        name=\$(basename "\${src}")
        cp -rL "\${src}" "best_models/\${name}"
        echo "  \${src} -> best_models/\${name}/"
    done

    # rmsd_effector_receptor_aligned is the whole-fit metric:
    # sequence-aligned single-pass Kabsch fit of the receptor, then that
    # transform applied rigidly to the effector.  Measures "once the
    # whole receptor is lined up, how far is the effector from where it
    # should be?" — the primary docking-quality question.  Lower is
    # better; null values sort to the end.
    #
    # Ties broken by ipsae_min descending (more confident first), then
    # model name and pdb_path for stable ordering.
    #
    # all_metrics.csv is kept in its original (model/seed) order so
    # downstream analysis scripts with positional assumptions still work.
    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python - << 'PYEOF'
import csv

in_path  = "all_metrics.csv"
out_path = "all_metrics_ranked_by_effector_rmsd.csv"

with open(in_path) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    rows = list(reader)

def sort_key(r):
    # Primary: whole-fit effector RMSD (lower is better).
    # Missing/empty → +inf so those rows sink to the bottom.
    try:
        eff = float(r.get("rmsd_effector_receptor_aligned") or "inf")
    except ValueError:
        eff = float("inf")
    # Secondary: ipsae_min descending (higher is better → negate).
    try:
        ipsae = -float(r.get("ipsae_min") or "0")
    except ValueError:
        ipsae = 0.0
    # Tertiary: model name + pdb_path for stable ordering.
    model = r.get("model", "")
    pdb   = r.get("pdb_path", "")
    return (eff, ipsae, model, pdb)

rows_sorted = sorted(rows, key=sort_key)

with open(out_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    writer.writeheader()
    for r in rows_sorted:
        writer.writerow(r)

print(f"Wrote {len(rows_sorted)} rows to {out_path}")
print(f"Top 10 by rmsd_effector_receptor_aligned (whole-fit):")
for r in rows_sorted[:10]:
    print(f"  {r.get('model', '?'):25s}  "
          f"eff_rmsd={r.get('rmsd_effector_receptor_aligned', ''):>6}  "
          f"ipsae_min={r.get('ipsae_min', ''):>6}  "
          f"{r.get('pdb_path', '')}")
PYEOF

    echo ""
    echo "=== Confidence Metrics Summary ==="
    column -t -s, all_metrics.csv | head -40 || cat all_metrics.csv
    """
}
