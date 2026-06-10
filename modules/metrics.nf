/*
 * =============================================================================
 * Metrics module — shared compute_metrics.py wrapper
 * =============================================================================
 * A single reusable process that runs compute_metrics.py against a
 * prediction directory and writes metrics.csv + best_model/.  Imported
 * with aliases in main.nf so each model gets its own results subdirectory:
 *
 *   include { COMPUTE_METRICS as METRICS_BOLTZ1 } from './modules/metrics'
 *   include { COMPUTE_METRICS as METRICS_BOLTZ2 } from './modules/metrics'
 *   ...
 *
 * Each alias is a separate process instance with its own work directory,
 * so they can run in parallel while still writing to distinct
 * ${outdir}/<model>/ publish targets.
 *
 * After compute_metrics.py writes best_model/, this process renames
 * best_model → ${model_tag}_best so the aggregation step can ingest a
 * flat list of uniquely-named directories with no tag/file bookkeeping.
 *
 * The model_tag value drives:
 *   - the --model argument passed to compute_metrics.py IF it matches a
 *     known parser, otherwise parser_tag is used
 *   - the publishDir subdirectory
 *   - the 'model' column in metrics.csv (rewritten post-hoc if parser_tag
 *     differs from model_tag, so variants like boltz1_msa are
 *     distinguishable in the aggregated CSV)
 */


/*
 * COMPUTE_METRICS
 * ---------------
 * Input tuple:
 *   model_tag      — identifier used for output subdir and the CSV model column
 *                    (e.g. boltz1, boltz1_msa, boltz2_constrained, af3_nomsa)
 *   parser_tag     — base parser to pass to --model (e.g. boltz1, boltz2,
 *                    chai1, af2m, af3, colabfold).  Normalises variants
 *                    back to their underlying output format.
 *   prediction_dir — staging dir containing all seed outputs for this model
 *   reference_pdb  — reference structure for RMSD/DockQ
 *   chain_a_len    — receptor chain length (int)
 *   chain_b_len    — effector chain length (int)
 *   rec_chain      — receptor chain ID (e.g. A, B)
 *   eff_chain      — effector chain ID (e.g. B, C)
 *
 * Outputs two tagged tuples so the aggregate step can correlate metrics
 * and best-model dirs by model_tag without extra plumbing.
 */
process COMPUTE_METRICS {
    tag "${model_tag}"
    label 'cpu'

    publishDir "${params.outdir}/${model_tag}", mode: 'copy',
        pattern: '{metrics.csv,best_model/**}'

    input:
    tuple val(model_tag),
          val(parser_tag),
          path(prediction_dir, stageAs: 'predictions'),
          path(reference_pdb),
          val(chain_a_len),
          val(chain_b_len),
          val(rec_chain),
          val(eff_chain)

    output:
    // For aggregation: renamed metrics CSV and renamed best-model dir.
    // The *_metrics.csv and *_best dir names are unique across all
    // model aliases, eliminating naming collisions when they get
    // staged side-by-side in AGGREGATE_RESULTS.
    path "${model_tag}_metrics.csv",  emit: tagged_metrics
    path "${model_tag}_best",         emit: tagged_best_dir
    // For per-model publishDir: the non-prefixed versions.
    path "metrics.csv",               emit: metrics
    path "best_model",                emit: best_model

    script:
    """
    set -euo pipefail

    mkdir -p best_model

    singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
        python ${projectDir}/bin/compute_metrics.py \\
            --model ${parser_tag} \\
            --prediction-dir predictions \\
            --chain-lengths ${chain_a_len} ${chain_b_len} \\
            --output-csv metrics.csv \\
            --best-model-dir best_model \\
            --reference-pdb ${reference_pdb} \\
            --receptor-chain ${rec_chain} \\
            --effector-chain ${eff_chain}

    # ─── Rewrite the 'model' column if parser_tag != model_tag ───────────
    # Needed so boltz1_msa / boltz1_constrained / af3_nomsa / colabfold_nomsa
    # appear as themselves in the aggregated CSV rather than as their base
    # parser name.
    if [ "${parser_tag}" != "${model_tag}" ]; then
        singularity exec --bind \${PWD}:\${PWD} ${params.benchmark_container} \\
            python - << 'PYEOF'
import csv
rows = []
with open("metrics.csv") as f:
    reader = csv.reader(f)
    header = next(reader)
    rows.append(header)
    try:
        mcol = header.index("model")
    except ValueError:
        mcol = 0
    for row in reader:
        if row:
            row[mcol] = "${model_tag}"
        rows.append(row)
with open("metrics.csv", "w", newline="") as f:
    csv.writer(f).writerows(rows)
PYEOF
    fi

    # ─── Create uniquely-named copies for the aggregate step ─────────────
    # These live alongside the non-prefixed outputs; publishDir only
    # copies the non-prefixed ones (via 'pattern').
    cp metrics.csv "${model_tag}_metrics.csv"
    cp -r best_model "${model_tag}_best"
    """
}
