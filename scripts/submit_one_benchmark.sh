#!/bin/bash
# Submit a SINGLE benchmark experiment by PDB id, from within its own directory
# (so SLURM logs and results land inside experiments/benchmarks/<PDB>/).
#
# Use this for a controlled test run instead of submit_all_benchmarks.sh, which
# fires all 43 Nextflow drivers at once and is hard to observe/debug.
#
#   Usage:  bash scripts/submit_one_benchmark.sh <PDB_ID>
#   e.g.    bash scripts/submit_one_benchmark.sh 5A6W
#
# Tip: for a faster smoke test, trim the `models:` list in the benchmark's
# params.yml to a couple of GPU models (e.g. boltz2, chai1) before submitting.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"

PDB="${1:-}"
if [ -z "${PDB}" ]; then
    echo "Usage: bash scripts/submit_one_benchmark.sh <PDB_ID>"
    echo "Available benchmarks:"
    ls "${PIPELINE_DIR}/experiments/benchmarks/" | tr '\n' ' '; echo
    exit 1
fi

DIR="${PIPELINE_DIR}/experiments/benchmarks/${PDB}"
if [ ! -d "${DIR}" ]; then
    echo "ERROR: no benchmark dir for '${PDB}' at ${DIR}"
    exit 1
fi
if [ ! -f "${DIR}/params.yml" ]; then
    echo "ERROR: ${DIR}/params.yml not found"
    exit 1
fi

JOBID="$(cd "${DIR}" && sbatch --parsable "${PIPELINE_DIR}/run_benchmark.slurm.sh" params.yml)"

cat <<EOF
Submitted ${PDB} as SLURM job ${JOBID}.

Watch the Nextflow driver:
  tail -f ${DIR}/.nextflow.log
  tail -f ${DIR}/nextflow_benchmark_${JOBID}.out

Watch what it submits to the GPU queue:
  watch -n 30 'squeue -u \$USER'

If a GPU child job ends oddly, get the verdict (State/Reason/Elapsed):
  sacct -j <child_jobid> --format=JobID,JobName%20,State,Elapsed,ExitCode,Reason
EOF
