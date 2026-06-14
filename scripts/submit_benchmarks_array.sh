#!/bin/bash
# Submit INCOMPLETE benchmarks as a SLURM job array, throttled to K concurrent,
# with auto-requeue on node failure. Replaces the controller approach: there is
# no babysitter process to die — SLURM enforces the %K concurrency and
# re-queues tasks if a node is rebooted/fails mid-run (e.g. GPU maintenance).
#
# Selection: choose benchmarks by tier (from the "Tier N" line in each
# params.yml), then drop any that are already complete (have a
# <PDB>_benchmark_results/all_metrics.csv) so a re-run only does outstanding work.
#
#   bash scripts/submit_benchmarks_array.sh --tiers 2,3,4          # default 6 concurrent
#   bash scripts/submit_benchmarks_array.sh --tiers 2,3,4 -j 4
#   bash scripts/submit_benchmarks_array.sh --tiers 2,3,4 --list   # preview, submit nothing
#   bash scripts/submit_benchmarks_array.sh --tiers all --include-complete  # don't skip completed
#
# The frozen PDB list and the array's SLURM logs go under experiments/array_results/,
# whose name matches sync_to_hpc.sh's *_results/ exclude (so a sync can't delete
# them mid-run).

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BENCH_DIR="${PIPELINE_DIR}/experiments/benchmarks"
RUNNER="${PIPELINE_DIR}/scripts/run_benchmarks_array.slurm.sh"
STATE_DIR="${PIPELINE_DIR}/experiments/array_results"

MAX_CONCURRENT=6
TIERS=""
LIST_ONLY=0
INCLUDE_COMPLETE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)          MAX_CONCURRENT="$2"; shift 2 ;;
        --tiers)            TIERS="$2"; shift 2 ;;
        --list|--dry)       LIST_ONLY=1; shift ;;
        --include-complete) INCLUDE_COMPLETE=1; shift ;;
        -h|--help)          sed -n '2,22p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ -z "${TIERS}" ]; then
    echo "ERROR: specify --tiers (e.g. --tiers 1 | --tiers 2,3,4 | --tiers all)" >&2
    exit 2
fi

if [ "${TIERS}" = "all" ]; then
    TIER_RE='Tier ?[1-9][: ]'
else
    digits="$(echo "${TIERS}" | tr -cd '0-9')"
    [ -n "${digits}" ] || { echo "ERROR: --tiers must be 'all' or digits" >&2; exit 2; }
    TIER_RE="Tier ?[${digits}][: ]"
fi

# Resolve tier set, then keep only the incomplete ones.
TODO=()
SKIP=()
while IFS= read -r pdb; do
    [ -n "${pdb}" ] || continue
    if [ "${INCLUDE_COMPLETE}" -eq 0 ] && \
       [ -f "${BENCH_DIR}/${pdb}/${pdb}_benchmark_results/all_metrics.csv" ]; then
        SKIP+=("${pdb}"); continue
    fi
    TODO+=("${pdb}")
done < <(grep -lE "${TIER_RE}" "${BENCH_DIR}"/*/params.yml 2>/dev/null \
            | xargs -n1 dirname | xargs -n1 basename | sort)

if [ "${#TODO[@]}" -eq 0 ]; then
    echo "Nothing to run for tiers '${TIERS}' (all matching benchmarks already complete)."
    [ "${#SKIP[@]}" -gt 0 ] && echo "Complete (skipped): ${SKIP[*]}"
    exit 0
fi

echo "To run (tiers ${TIERS}, ${#TODO[@]}): ${TODO[*]}"
[ "${#SKIP[@]}" -gt 0 ] && echo "Already complete, skipped (${#SKIP[@]}): ${SKIP[*]}"
echo "Max concurrent: ${MAX_CONCURRENT}   (auto-requeue on node failure)"

if [ "${LIST_ONLY}" -eq 1 ]; then
    echo "(--list: nothing submitted)"
    exit 0
fi

mkdir -p "${STATE_DIR}"
stamp="$(date +%Y%m%d_%H%M%S)"
LISTFILE="${STATE_DIR}/pdb_list_${stamp}.txt"
printf '%s\n' "${TODO[@]}" > "${LISTFILE}"
N="${#TODO[@]}"

jobid="$(sbatch --parsable \
    --array="0-$((N - 1))%${MAX_CONCURRENT}" \
    --output="${STATE_DIR}/nf_%A_%a.out" \
    --error="${STATE_DIR}/nf_%A_%a.err" \
    "${RUNNER}" "${LISTFILE}")"

echo "Submitted array job ${jobid}: ${N} tasks, %${MAX_CONCURRENT} concurrent, auto-requeue."
echo "List file: ${LISTFILE}"
echo "Watch:  squeue -u \$USER"
echo "Logs:   ${STATE_DIR}/nf_${jobid}_*.out"
