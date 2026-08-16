#!/bin/bash
#SBATCH --job-name="spb_pytest"
#SBATCH -p jic-short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=pytest_%j.out
#SBATCH --error=pytest_%j.err

# Run the test suite on the HPC.
#
# Three tiers, one marker each:
#   local_unit         pure helpers on synthetic input, runs anywhere
#   local_integration  reads the committed references and metrics
#   hpc                needs a GPU, the containers, or a finished run
#
# The hpc tier reads published best_models/ trees, so run it after a benchmark.
# Point it elsewhere with SPB_BENCHMARKS_DIR. Tests skip cleanly when their
# inputs are absent rather than failing.
#
# Default marker is the full suite; override with MARKER, for example
#   MARKER=hpc sbatch tests/run_pytest.slurm.sh
#
# Usage:
#   sbatch tests/run_pytest.slurm.sh                 # full suite (yaml tests skip)
#   sbatch tests/run_pytest.slurm.sh -m local_unit   # one tier
#   sbatch tests/run_pytest.slurm.sh tests/unit -q   # any pytest args

set -euo pipefail

REPO_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"
CONTAINER="${PYTEST_CONTAINER:-/hpc-home/jowillia/singularity/pytest/pytest_runner.img}"

cd "${REPO_DIR}"

# Default to the whole suite if no pytest args are supplied.
if [ "$#" -eq 0 ]; then
    set -- tests/ -ra
fi

echo "============================================================"
echo "structure-prediction-benchmarking — pytest"
echo "Repo:      ${REPO_DIR}"
echo "Container: ${CONTAINER}"
echo "Args:      $*"
echo "Date:      $(date)"
echo "Node:      $(hostname)"
echo "============================================================"

singularity exec --bind "${REPO_DIR}:${REPO_DIR}" "${CONTAINER}" \
    python -m pytest "$@"
