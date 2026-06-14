#!/bin/bash
# Run the throttled submitter for the OTHER 25 benchmarks (Tiers 2, 3 and 4) AS
# A SLURM JOB, so you don't have to stay logged in.  Same approach as
# submit_tier1_controller.slurm.sh — this tiny "controller" job sits on a compute
# node, submits the drivers in a rolling window (default 6 concurrent), then
# exits.  Survives logout; nothing runs on the login node.
#
#   Usage:
#     sbatch scripts/submit_rest_controller.slurm.sh           # default 6 concurrent
#     CONC=4 sbatch scripts/submit_rest_controller.slurm.sh    # cap at 4 concurrent
#
#   Watch:
#     squeue -u $USER                       # rest_ctl = this controller; nf_benchmark = drivers
#     cat rest_controller_<jobid>.out       # controller progress log
#
#SBATCH --job-name=rest_ctl
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=512M
#SBATCH --time=48:00:00
#SBATCH --output=rest_controller_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

set -euo pipefail

PIPELINE_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"
CONC="${CONC:-6}"

echo "============================================================"
echo "Tiers 2/3/4 controller starting on $(hostname) at $(date)"
echo "Max concurrent drivers: ${CONC}"
echo "============================================================"

# Rolling-window submission of the other 25 (Tiers 2,3,4). Runs here on a
# compute node, not the login node, so it keeps going after you log out.
bash "${PIPELINE_DIR}/scripts/submit_benchmarks.sh" --tiers 2,3,4 -j "${CONC}"

echo "============================================================"
echo "Controller done at $(date) — all Tier 2/3/4 drivers submitted."
echo "Drivers will continue/finish on their own; this job now exits."
echo "============================================================"
