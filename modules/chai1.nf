/*
 * =============================================================================
 * Chai-1 module
 * =============================================================================
 * Chai-1 tries to write weights to CHAI_DOWNLOADS_DIR at runtime.  The
 * container already ships pre-downloaded weights at /opt/chai_cache.
 * The strategy (preserved from bench_chai1.sh):
 *
 *   1. Create a small job-local scratch dir under /tmp for Chai-1
 *      metadata writes (not the multi-GB weights).
 *   2. Expose /opt/chai_cache read-only via CHAI_DOWNLOADS_DIR so Chai-1
 *      finds weights immediately without downloading anything.
 *   3. Pre-populate the scratch with tiny top-level metadata files from
 *      /opt/chai_cache so Chai-1 doesn't attempt a network download for
 *      index files (the HPC has no internet access).
 *   4. Clean up the scratch dir on EXIT via a trap.
 *
 * Seeds are run SERIALLY in a single Python heredoc because Chai-1's
 * run_inference() is a Python function rather than a CLI.
 */


/*
 * CHAI1
 * -----
 * Run Chai-1 v0.6.1 inference across 5 seeds × 5 diffusion samples.
 */
process CHAI1 {
    tag "${params.project_name}"
    label 'gpu'

    publishDir "${params.outdir}/chai1", mode: 'copy',
        pattern: '{input.fasta,all_outputs/**}'

    input:
    val receptor_seq
    val effector_seq
    val rec_chain
    val eff_chain

    output:
    path "all_outputs",  emit: prediction_dir
    path "input.fasta"

    script:
    """
    set -euo pipefail

    # ─── Chai-1 FASTA (uses protein|name=<header> style) ─────────────────
    cat > input.fasta << EOF
>protein|name=chain_A
${receptor_seq}
>protein|name=chain_B
${effector_seq}
EOF

    # ─── Job-local scratch cache for Chai-1 metadata writes ───────────────
    CHAI_SCRATCH=\$(mktemp -d "/tmp/chai1_cache_\${SLURM_JOB_ID:-\$\$}_XXXXXX")
    echo "Chai-1 scratch cache: \${CHAI_SCRATCH}"

    cleanup_chai_cache() {
        if [ -d "\${CHAI_SCRATCH}" ]; then
            echo "Removing Chai-1 scratch cache (\${CHAI_SCRATCH})..."
            rm -rf "\${CHAI_SCRATCH}"
            echo "  Done."
        fi
    }
    trap cleanup_chai_cache EXIT

    # Pre-populate scratch with top-level metadata from the container cache
    # so Chai-1 doesn't attempt a network download for index files.
    singularity exec ${params.benchmark_container} \\
        bash -c 'find /opt/chai_cache -maxdepth 1 -type f 2>/dev/null | head -20' \\
        | while read -r f; do
            singularity exec ${params.benchmark_container} cat "\${f}" > "\${CHAI_SCRATCH}/\$(basename "\${f}")" 2>/dev/null || true
        done

    echo "Running Chai-1 complex prediction..."

    export CHAI_WORKDIR="\${PWD}"

    singularity exec --nv \\
        --bind \${PWD}:\${PWD} \\
        --bind "\${CHAI_SCRATCH}:/chai_scratch:rw" \\
        --env CHAI_DOWNLOADS_DIR="/opt/chai_cache" \\
        --env CHAI_WORKDIR="\${PWD}" \\
        ${params.benchmark_container} python << 'CHAI_SCRIPT'
import os
from pathlib import Path

workdir = os.environ["CHAI_WORKDIR"]
output_dir = Path(workdir) / "output"
output_dir.mkdir(parents=True, exist_ok=True)
fasta_path = Path(workdir) / "input.fasta"

print(f"Input:  {fasta_path}")
print(f"Output: {output_dir}")
print(f"CHAI_DOWNLOADS_DIR: {os.environ.get('CHAI_DOWNLOADS_DIR', 'not set')}")

from chai_lab.chai1 import run_inference

result = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    num_trunk_recycles=20,
    num_diffn_timesteps=20,
    num_diffn_samples=5,
    seed=42,
)

print("\\n--- Chai-1 Results (seed 42) ---")
for attr in ['aggregate_score', 'ptm', 'iptm', 'ranking_confidence', 'plddt', 'scores']:
    if hasattr(result, attr):
        print(f"  {attr} = {getattr(result, attr)}")

# Run 4 more seeds serially
for seed in [123, 456, 789, 1024]:
    seed_dir = Path(workdir) / f"output_seed{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)
    print(f"\\n--- Running Chai-1 with seed {seed} ---")
    result = run_inference(
        fasta_file=fasta_path,
        output_dir=seed_dir,
        num_trunk_recycles=20,
        num_diffn_timesteps=20,
        num_diffn_samples=5,
        seed=seed,
    )
    print(f"--- Chai-1 Results (seed {seed}) ---")
    for attr in ['aggregate_score', 'ptm', 'iptm', 'ranking_confidence', 'plddt', 'scores']:
        if hasattr(result, attr):
            print(f"  {attr} = {getattr(result, attr)}")

print("\\nOutput files:")
for f in sorted(output_dir.rglob("*")):
    if f.is_file():
        print(f"  {f} ({f.stat().st_size} bytes)")
CHAI_SCRIPT

    # ─── Aggregate multi-seed outputs ─────────────────────────────────────
    bash ${projectDir}/bin/aggregate_seed_outputs.sh cif pdb npz pt json

    echo "Aggregated PDB files:"
    find all_outputs -name "*.pdb" | sort | head -20 || true
    """
}
