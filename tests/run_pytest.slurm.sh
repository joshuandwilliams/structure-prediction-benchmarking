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

# Run the Python test suite on the HPC, inside the lightweight pytest_runner
# container (pytest + numpy + pandas). That covers the whole local_unit tier:
# the numpy maths, the pure-Python helpers, and the extractor/FASTA CLIs (which
# are pure stdlib). pytest_runner does NOT ship gemmi/biopython/PyYAML, so the
# PyYAML-dependent validator tests skip cleanly (they importorskip "yaml").
#
# A genuine `hpc`-tier run that parses real predictions with gemmi would need
# the benchmark container instead — override CONTAINER below for that.
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
