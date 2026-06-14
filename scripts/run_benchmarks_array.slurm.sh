#!/bin/bash
# Job-array task runner: one array index = one benchmark.
#
# Submitted by submit_benchmarks_array.sh, which sets --array=0-N%K and passes a
# frozen list file of PDB ids. This task maps its SLURM_ARRAY_TASK_ID to a PDB,
# cds into that benchmark dir, and runs the normal pipeline launcher.
#
# --requeue: if the node is rebooted / fails (e.g. GPU maintenance), SLURM
# automatically re-queues this task instead of failing it. The launcher runs
# Nextflow with -resume, so a re-queued task continues from its last cached
# step rather than restarting from scratch. No controller/babysitter process is
# involved — SLURM itself manages throttling (%K) and re-queueing.
#
#SBATCH --job-name=nf_benchmark
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --requeue
#SBATCH --output=array_%A_%a.out   # overridden by submit_benchmarks_array.sh

set -euo pipefail

PIPELINE_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"

LISTFILE="${1:?usage: sbatch --array=0-N%K run_benchmarks_array.slurm.sh <pdb_list_file>}"
idx="${SLURM_ARRAY_TASK_ID:?must be run as a SLURM job array}"

PDB="$(sed -n "$((idx + 1))p" "$LISTFILE")"
if [ -z "${PDB}" ]; then
    echo "ERROR: no PDB on line $((idx + 1)) of ${LISTFILE}" >&2
    exit 1
fi

DIR="${PIPELINE_DIR}/experiments/benchmarks/${PDB}"
if [ ! -f "${DIR}/params.yml" ]; then
    echo "ERROR: ${DIR}/params.yml not found" >&2
    exit 1
fi

echo "============================================================"
echo "Array task ${idx} -> ${PDB}"
echo "Restart count: ${SLURM_RESTART_COUNT:-0}   Node: $(hostname)   $(date)"
echo "============================================================"

cd "${DIR}"
exec bash "${PIPELINE_DIR}/run_benchmark.slurm.sh" params.yml
