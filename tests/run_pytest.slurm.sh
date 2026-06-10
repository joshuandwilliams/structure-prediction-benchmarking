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

# Run the Python test suite on the HPC, inside the benchmark container so the
# parser dependencies (gemmi, numpy, biopython, PyYAML) are available.
#
# Usage:
#   sbatch tests/run_pytest.slurm.sh                 # full suite
#   sbatch tests/run_pytest.slurm.sh -m local_unit   # one tier
#   sbatch tests/run_pytest.slurm.sh tests/unit -q   # any pytest args

set -euo pipefail

REPO_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"
BENCHMARK_IMG="/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/Boltz1_Boltz2_Chai1.img"

cd "${REPO_DIR}"

# Default to the whole suite if no pytest args are supplied.
if [ "$#" -eq 0 ]; then
    set -- tests/ -ra
fi

echo "============================================================"
echo "structure-prediction-benchmarking — pytest"
echo "Repo:      ${REPO_DIR}"
echo "Container: ${BENCHMARK_IMG}"
echo "Args:      $*"
echo "Date:      $(date)"
echo "Node:      $(hostname)"
echo "============================================================"

singularity exec --bind "${REPO_DIR}:${REPO_DIR}" "${BENCHMARK_IMG}" \
    python -m pytest "$@"
