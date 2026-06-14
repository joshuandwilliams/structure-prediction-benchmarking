#!/bin/bash
# Submit a selected set of benchmarks, throttled to a fixed number of concurrent
# Nextflow drivers (rolling window), so we never overload `squeue` the way 43
# simultaneous drivers did.
#
# Generalised version of submit_tier1_benchmarks.sh: choose which benchmarks by
# tier, read from the "Tier N" annotation in each params.yml.
#
#   --tiers 1          # just Tier 1 (the 18 HMA targets)
#   --tiers 2,3,4      # Tiers 2, 3 and 4 (the other 25)
#   --tiers all        # everything (43)
#
# Each benchmark is one head job (run_benchmark.slurm.sh, job name
# "nf_benchmark") submitted from inside its own experiments/benchmarks/<PDB>/
# directory. At most MAX_CONCURRENT (default 6) drivers run at once; the next is
# submitted as soon as a slot frees.
#
#   Usage:
#     bash scripts/submit_benchmarks.sh --tiers 2,3,4          # default 6 concurrent
#     bash scripts/submit_benchmarks.sh --tiers 2,3,4 -j 4     # cap at 4
#     bash scripts/submit_benchmarks.sh --tiers 2,3,4 --list   # preview, submit nothing
#
# It sleeps while waiting for slots, so run it under tmux/nohup — or use the
# matching *_controller.slurm.sh wrapper to run it as a SLURM job (hands-off).

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BENCH_DIR="${PIPELINE_DIR}/experiments/benchmarks"
JOB_NAME="nf_benchmark"          # head-job name set in run_benchmark.slurm.sh

MAX_CONCURRENT=6
POLL_INTERVAL=60                 # seconds between slot checks
LIST_ONLY=0
TIERS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)    MAX_CONCURRENT="$2"; shift 2 ;;
        --tiers)      TIERS="$2"; shift 2 ;;
        --list|--dry) LIST_ONLY=1; shift ;;
        -h|--help)    sed -n '2,30p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ -z "${TIERS}" ]; then
    echo "ERROR: specify --tiers (e.g. --tiers 1 | --tiers 2,3,4 | --tiers all)" >&2
    exit 2
fi

# ── Build a tier-matching regex from the selection ────────────────────────
if [ "${TIERS}" = "all" ]; then
    TIER_RE='Tier ?[1-9][: ]'
else
    digits="$(echo "${TIERS}" | tr -cd '0-9')"   # "2,3,4" -> "234"
    if [ -z "${digits}" ]; then
        echo "ERROR: --tiers must be 'all' or digits like 1 or 2,3,4" >&2
        exit 2
    fi
    TIER_RE="Tier ?[${digits}][: ]"
fi

# ── Resolve the selected set (sorted, deterministic) ──────────────────────
SET=()
while IFS= read -r pdb; do
    [ -n "${pdb}" ] && SET+=("${pdb}")
done < <(
    grep -lE "${TIER_RE}" "${BENCH_DIR}"/*/params.yml 2>/dev/null \
        | xargs -n1 dirname | xargs -n1 basename | sort
)

if [ "${#SET[@]}" -eq 0 ]; then
    echo "ERROR: no benchmarks matched tiers '${TIERS}' under ${BENCH_DIR}" >&2
    exit 1
fi

echo "Selected (tiers ${TIERS}, ${#SET[@]}): ${SET[*]}"
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
for pdb in "${SET[@]}"; do
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
echo "All ${#SET[@]} benchmarks (tiers ${TIERS}) submitted."
echo "Watch drivers:  squeue -u \$USER -n ${JOB_NAME}"
echo "Watch a run:    tail -f ${BENCH_DIR}/<PDB>/.nextflow.log"
