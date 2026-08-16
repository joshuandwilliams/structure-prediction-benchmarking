#!/bin/bash
# Run a set of benchmarks as a throttled, auto-requeueing SLURM job array.
#
# This script is BOTH the launcher and the array task. It branches on whether
# SLURM_ARRAY_TASK_ID is set:
#
#   no  -> launcher mode. Runs on the login node. Resolves which targets to
#          run, freezes the list to a file, and re-submits ITSELF with
#          `sbatch --array=0-N%K`. The array size N is only knowable after the
#          scan, which is why a launcher pass has to exist at all.
#   yes -> array-task mode. One task per target: map the index to a PDB, cd
#          into its directory, and hand off to the Nextflow driver.
#
# Selection: choose benchmarks by tier (from the "Tier N" line in the chosen
# params file), then drop any already complete. Completion is the
# all_metrics.csv under that params file's own `outdir:`, so Run 1 and Run 2
# are tracked independently.
#
# --params picks which of a target's params files to run. Each Tier-1 target
# carries two: params.yml (Run 1, the 13 template-free variants) and
# params_template.yml (Run 2, the three effector-template Boltz-2 variants).
#
#   bash scripts/run_benchmarks.sh --tiers 1                              # Run 1
#   bash scripts/run_benchmarks.sh --tiers 1 --params params_template.yml # Run 2
#   bash scripts/run_benchmarks.sh --tiers 1 -j 4
#   bash scripts/run_benchmarks.sh --tiers 1 --pdbs 6G10   # one target only
#   bash scripts/run_benchmarks.sh --tiers 1 --params params_template.yml \
#       --after 12345                                      # chain Run 2 after Run 1
#   bash scripts/run_benchmarks.sh --tiers 1 --list    # preview, submit nothing
#   bash scripts/run_benchmarks.sh --tiers 1 --include-complete
#
# Submit Run 2 only after Run 1, and from the same target directories: the
# Nextflow cache then satisfies the shared COLABFOLD_SEARCH step from Run 1
# rather than recomputing it (~34 min per target).
#
# The driver is invoked with `exec bash`, NOT sbatch: the array task IS the
# allocation the Nextflow driver should run inside. Submitting a nested job
# would make the task exit in seconds, so SLURM's %K would throttle
# submissions rather than concurrent runs, and --requeue would have nothing
# left to requeue. run_benchmark.slurm.sh's own #SBATCH headers are inert
# under bash; the headers below are the ones that take effect, and they match.
#
# --requeue: if a node is rebooted or fails mid-run (e.g. GPU maintenance),
# SLURM re-queues the task instead of failing it. The driver runs Nextflow
# with -resume, so a re-queued task continues from its last cached step. No
# babysitter process is involved — SLURM manages throttling and re-queueing.
#
# The frozen PDB list and the array's SLURM logs go under
# experiments/array_results/, whose name matches sync_to_hpc.sh's *_results/
# exclude, so a sync cannot delete them mid-run.
#
#SBATCH --job-name=nf_benchmark
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --requeue
#SBATCH --output=array_%A_%a.out

set -euo pipefail

SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
PIPELINE_DIR="$(cd "$(dirname "${SELF}")/.." && pwd)"
BENCH_DIR="${PIPELINE_DIR}/experiments/benchmarks"
STATE_DIR="${PIPELINE_DIR}/experiments/array_results"


# ═════════════════════════════════════════════════════════════════════════
# Array-task mode
# ═════════════════════════════════════════════════════════════════════════
# Positional args come from the launcher's sbatch call below.

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    LISTFILE="${1:?array mode: missing pdb list file}"
    PARAMS="${2:-params.yml}"
    # SLURM copies the submitted script into its own spool directory, so $0
    # points at /var/spool/slurmd/... rather than into the repo and the
    # PIPELINE_DIR derived above is wrong. The launcher knows the real root and
    # passes it as $3; SLURM_SUBMIT_DIR is the fallback for a hand-built sbatch.
    PIPELINE_DIR="${3:-${SLURM_SUBMIT_DIR:-${PIPELINE_DIR}}}"
    BENCH_DIR="${PIPELINE_DIR}/experiments/benchmarks"
    idx="${SLURM_ARRAY_TASK_ID}"

    if [ ! -d "${BENCH_DIR}" ]; then
        echo "ERROR: ${BENCH_DIR} is not a directory. Array mode needs the" >&2
        echo "       repo root as \$3, or SLURM_SUBMIT_DIR set." >&2
        exit 1
    fi

    PDB="$(sed -n "$((idx + 1))p" "${LISTFILE}")"
    if [ -z "${PDB}" ]; then
        echo "ERROR: no PDB on line $((idx + 1)) of ${LISTFILE}" >&2
        exit 1
    fi

    DIR="${BENCH_DIR}/${PDB}"
    if [ ! -f "${DIR}/${PARAMS}" ]; then
        echo "ERROR: ${DIR}/${PARAMS} not found" >&2
        exit 1
    fi

    echo "============================================================"
    echo "Array task ${idx} -> ${PDB}  (${PARAMS})"
    echo "Restart count: ${SLURM_RESTART_COUNT:-0}   Node: $(hostname)   $(date)"
    echo "============================================================"

    cd "${DIR}"
    exec bash "${PIPELINE_DIR}/run_benchmark.slurm.sh" "${PARAMS}"
fi


# ═════════════════════════════════════════════════════════════════════════
# Launcher mode
# ═════════════════════════════════════════════════════════════════════════

MAX_CONCURRENT=6
TIERS=""
PARAMS="params.yml"
LIST_ONLY=0
INCLUDE_COMPLETE=0
ONLY_PDBS=""
AFTER_JOB=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -j|--jobs)          MAX_CONCURRENT="$2"; shift 2 ;;
        --tiers)            TIERS="$2"; shift 2 ;;
        --params)           PARAMS="$2"; shift 2 ;;
        --pdbs)             ONLY_PDBS="$2"; shift 2 ;;
        --after)            AFTER_JOB="$2"; shift 2 ;;
        --list|--dry)       LIST_ONLY=1; shift ;;
        --include-complete) INCLUDE_COMPLETE=1; shift ;;
        -h|--help)          sed -n '2,50p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [ -z "${TIERS}" ]; then
    echo "ERROR: specify --tiers (the benchmark set is tier 1: --tiers 1)" >&2
    exit 2
fi

if [ "${TIERS}" = "all" ]; then
    TIER_RE='Tier ?[1-9][: ]'
else
    digits="$(echo "${TIERS}" | tr -cd '0-9')"
    [ -n "${digits}" ] || { echo "ERROR: --tiers must be 'all' or digits" >&2; exit 2; }
    TIER_RE="Tier ?[${digits}][: ]"
fi

# The results directory this params file writes to, read from its own
# `outdir:` line, so the completion check follows the run rather than assuming
# the Run-1 naming.
outdir_of() {
    sed -n 's/^outdir:[[:space:]]*"\{0,1\}\.\{0,1\}\/\{0,1\}\([^"]*\)"\{0,1\}[[:space:]]*$/\1/p' \
        "$1" | head -1
}

# Resolve tier set, then keep only the incomplete ones.
TODO=()
SKIP=()
while IFS= read -r pdb; do
    [ -n "${pdb}" ] || continue
    if [ "${INCLUDE_COMPLETE}" -eq 0 ]; then
        out="$(outdir_of "${BENCH_DIR}/${pdb}/${PARAMS}")"
        if [ -n "${out}" ] && [ -f "${BENCH_DIR}/${pdb}/${out}/all_metrics.csv" ]; then
            SKIP+=("${pdb}"); continue
        fi
    fi
    TODO+=("${pdb}")
done < <(grep -lE "${TIER_RE}" "${BENCH_DIR}"/*/"${PARAMS}" 2>/dev/null \
            | xargs -n1 dirname | xargs -n1 basename | sort)

# --pdbs narrows the tier set to named targets, for a smoke test before
# committing the whole array. Validated against the scanned set so a typo
# fails here rather than silently running nothing.
if [ -n "${ONLY_PDBS}" ]; then
    requested="$(echo "${ONLY_PDBS}" | tr ',' ' ')"
    FILTERED=()
    for want in ${requested}; do
        found=0
        for have in "${TODO[@]+"${TODO[@]}"}" "${SKIP[@]+"${SKIP[@]}"}"; do
            [ "${want}" = "${have}" ] && { found=1; break; }
        done
        [ "${found}" -eq 1 ] || {
            echo "ERROR: ${want} is not a tier-${TIERS} target with ${PARAMS}" >&2
            exit 2
        }
        for have in "${TODO[@]+"${TODO[@]}"}"; do
            [ "${want}" = "${have}" ] && FILTERED+=("${want}")
        done
    done
    TODO=("${FILTERED[@]+"${FILTERED[@]}"}")
fi

if [ "${#TODO[@]}" -eq 0 ]; then
    echo "Nothing to run for tiers '${TIERS}' with ${PARAMS} (all complete)."
    [ "${#SKIP[@]}" -gt 0 ] && echo "Complete (skipped): ${SKIP[*]}"
    exit 0
fi

echo "Params file: ${PARAMS}"
echo "To run (tiers ${TIERS}, ${#TODO[@]}): ${TODO[*]}"
[ "${#SKIP[@]}" -gt 0 ] && echo "Already complete, skipped (${#SKIP[@]}): ${SKIP[*]}"
echo "Max concurrent: ${MAX_CONCURRENT}   (auto-requeue on node failure)"

if [ "${LIST_ONLY}" -eq 1 ]; then
    echo "(--list: nothing submitted)"
    exit 0
fi

mkdir -p "${STATE_DIR}"
stamp="$(date +%Y%m%d_%H%M%S)"
LISTFILE="${STATE_DIR}/pdb_list_${PARAMS%.yml}_${stamp}.txt"
printf '%s\n' "${TODO[@]}" > "${LISTFILE}"
N="${#TODO[@]}"

# Re-submit this same file as the array task. SLURM_ARRAY_TASK_ID will be set
# in those jobs, so the branch at the top takes the array path.
# Run 2 launches Nextflow in the SAME target directories as Run 1, so the two
# cannot overlap: they would share work/ and .nextflow and corrupt each other.
# --after chains this submission behind an earlier job id. afterany rather than
# afterok, because a target that fails in Run 1 should not block Run 2 on the
# other seventeen.
DEP=()
if [ -n "${AFTER_JOB}" ]; then
    DEP=(--dependency="afterany:${AFTER_JOB}")
    echo "Chained: starts after job ${AFTER_JOB} finishes."
fi

jobid="$(sbatch --parsable \
    "${DEP[@]+"${DEP[@]}"}" \
    --array="0-$((N - 1))%${MAX_CONCURRENT}" \
    --output="${STATE_DIR}/nf_%A_%a.out" \
    --error="${STATE_DIR}/nf_%A_%a.err" \
    "${SELF}" "${LISTFILE}" "${PARAMS}" "${PIPELINE_DIR}")"

echo "Submitted array job ${jobid}: ${N} tasks, %${MAX_CONCURRENT} concurrent, auto-requeue."
echo "List file: ${LISTFILE}"
echo "Watch:  squeue -u \$USER"
echo "Logs:   ${STATE_DIR}/nf_${jobid}_*.out"
