#!/bin/bash
# Run the throttled Tier 1 submitter AS A SLURM JOB, so you don't have to stay
# logged in.  This tiny "controller" job sits on a compute node and submits the
# Tier 1 driver jobs in a rolling window (default 6 concurrent), continuing
# automatically as slots free, then exits.  Survives logout; nothing runs on the
# login node.
#
#   Usage:
#     sbatch scripts/submit_tier1_controller.slurm.sh           # default 6 concurrent
#     CONC=4 sbatch scripts/submit_tier1_controller.slurm.sh    # cap at 4 concurrent
#
#   Watch:
#     squeue -u $USER                       # tier1_ctl = this controller; nf_benchmark = drivers
#     cat tier1_controller_<jobid>.out      # controller progress log
#
#SBATCH --job-name=tier1_ctl
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=512M
#SBATCH --time=48:00:00
#SBATCH --output=tier1_controller_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

set -euo pipefail

PIPELINE_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"
CONC="${CONC:-6}"

echo "============================================================"
echo "Tier 1 controller starting on $(hostname) at $(date)"
echo "Max concurrent drivers: ${CONC}"
echo "============================================================"

# Hands off the rolling-window submission to the tested submitter. It polls
# squeue and only keeps ${CONC} nf_benchmark drivers active at once, so we never
# overload squeue. It runs here on a compute node, not the login node, so it
# keeps going after you log out.
bash "${PIPELINE_DIR}/scripts/submit_tier1_benchmarks.sh" -j "${CONC}"

echo "============================================================"
echo "Controller done at $(date) — all Tier 1 drivers submitted."
echo "Drivers will continue/finish on their own; this job now exits."
echo "============================================================"
