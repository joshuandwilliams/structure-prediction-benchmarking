#!/bin/bash
# Submit all 43 benchmark experiments, each from within its own directory
# so that SLURM logs and results land inside experiments/benchmarks/<PDB>/.

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

for d in "${PIPELINE_DIR}/experiments/benchmarks/"/*/; do
    (cd "$d" && sbatch "${PIPELINE_DIR}/run_benchmark.slurm.sh" params.yml)
done
