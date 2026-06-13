#!/bin/bash
# Submit the Tier 1 (HMA receptor / effector) benchmarks, throttled to a fixed
# number of concurrent Nextflow drivers.
#
# Why throttle: submitting all 43 drivers at once (submit_all_benchmarks.sh)
# had ~43 Nextflow head jobs polling `squeue` simultaneously, which triggered
# transient squeue dropouts — Nextflow then thought running jobs had vanished
# and stalled. Keeping only a handful of drivers active at once avoids that.
#
# Each benchmark is one head job (run_benchmark.slurm.sh, job name
# "nf_benchmark") submitted from inside its own experiments/benchmarks/<PDB>/
# directory. We keep at most MAX_CONCURRENT head jobs active, waiting for a
# slot to free before submitting the next.
#
#   Usage:
#     bash scripts/submit_tier1_benchmarks.sh           # submit (default 6 concurrent)
#     bash scripts/submit_tier1_benchmarks.sh -j 4      # cap at 4 concurrent drivers
#     bash scripts/submit_tier1_benchmarks.sh --list    # show the Tier 1 set, submit nothing
#
# The Tier 1 set is derived from the "Tier 1" annotation in each params.yml,
# so it stays correct if benchmarks are added/retiered.
#
# Tip: it sleeps while waiting for slots, so run it under tmux/nohup if your
# login session might drop.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BENCH_DIR="${PIPELINE_DIR}/experiments/benchmarks"
JOB_NAME="nf_benchmark"          # head-job name set in run_benchmark.slurm.sh

MAX_CONCURRENT=6
POLL_INTERVAL=60                 # seconds between slot checks
LIST_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)    MAX_CONCURRENT="$2"; shift 2 ;;
        --list|--dry) LIST_ONLY=1; shift ;;
        -h|--help)    sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

# ── Resolve the Tier 1 set from the params.yml headers ────────────────────
# (while-read loop instead of mapfile for bash 3.2 portability)
TIER1=()
while IFS= read -r pdb; do
    [ -n "${pdb}" ] && TIER1+=("${pdb}")
done < <(
    grep -lE 'Tier ?1[: ]' "${BENCH_DIR}"/*/params.yml 2>/dev/null \
        | xargs -n1 dirname | xargs -n1 basename | sort
)

if [ "${#TIER1[@]}" -eq 0 ]; then
    echo "ERROR: no Tier 1 benchmarks found under ${BENCH_DIR}" >&2
    exit 1
fi

echo "Tier 1 benchmarks (${#TIER1[@]}): ${TIER1[*]}"
echo "Max concurrent drivers: ${MAX_CONCURRENT}"

if [ "${LIST_ONLY}" -eq 1 ]; then
    echo "(--list: nothing submitted)"
    exit 0
fi

# ── Count this user's currently-active head jobs (R/PD/CG all count) ───────
active_heads() {
    squeue -u "$USER" -h -n "${JOB_NAME}" 2>/dev/null | wc -l | tr -d ' '
}

# ── Submit with a rolling concurrency cap ─────────────────────────────────
for pdb in "${TIER1[@]}"; do
    dir="${BENCH_DIR}/${pdb}"
    if [ ! -f "${dir}/params.yml" ]; then
        echo "[${pdb}] SKIP: no params.yml"
        continue
    fi

    while [ "$(active_heads)" -ge "${MAX_CONCURRENT}" ]; do
        echo "  $(date '+%H:%M:%S')  at cap ($(active_heads)/${MAX_CONCURRENT} ${JOB_NAME} active); waiting ${POLL_INTERVAL}s…"
        sleep "${POLL_INTERVAL}"
    done

    jobid="$(cd "${dir}" && sbatch --parsable "${PIPELINE_DIR}/run_benchmark.slurm.sh" params.yml)"
    echo "[${pdb}] submitted job ${jobid}  (active now: $(active_heads)/${MAX_CONCURRENT})"
done

echo
echo "All ${#TIER1[@]} Tier 1 benchmarks submitted."
echo "Watch drivers:  squeue -u \$USER -n ${JOB_NAME}"
echo "Watch a run:    tail -f ${BENCH_DIR}/<PDB>/.nextflow.log"
