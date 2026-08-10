#!/bin/bash
#
# Pull benchmark run outputs back from the HPC to the Mac for local analysis.
#
# The pipeline runs on the HPC and writes its results under
# experiments/<target>/<project>_results/.  Those trees are gitignored and
# never pushed up by sync_to_hpc.sh; this script pulls a chosen target's
# results down so they can be analysed locally (plotting, comparison) without
# committing the large prediction artefacts.
#
# Usage:
#   ./scripts/sync_from_hpc.sh --target 6G10            # one target, real sync
#   ./scripts/sync_from_hpc.sh --target 6G10 --dry      # dry-run
#   ./scripts/sync_from_hpc.sh --all                    # every target on HPC
#
# Transport:
#   rsync over SSH using the 'slurm' host alias defined in ~/.ssh/config.
#
# Scope:
#   Pulls experiments/<target>/ result trees only.  Source-tree code stays on
#   the Mac as authoritative.  Nextflow scratch (work/, .nextflow*) is
#   excluded, so only the published *_results/ outputs come back.
#
#   SLURM job logs (*.out, *.err) ARE pulled deliberately.  They record what
#   actually happened in each run and are the first thing to read when a
#   benchmark fails.  They are gitignored, so they stay local and untracked.
#
# Safety:
#   This script never passes --delete.  It can only add or update files on the
#   Mac, never remove them.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HPC_USER_HOST="slurm"
HPC_BASE="receptor_design/structure-prediction-benchmarking"

DRY_RUN=""
TARGET=""
ALL=""

usage() {
    cat <<EOF
Usage: $(basename "$0") (--target NAME | --all) [--dry]

  --target NAME  Pull experiments/NAME/ results from the HPC.
  --all          Pull every experiments/<target>/ that exists on the HPC.
  --dry          rsync dry-run; show what would change without copying.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry|--dry-run) DRY_RUN="--dry-run"; shift ;;
        --target)        TARGET="$2"; shift 2 ;;
        --all)           ALL="1"; shift ;;
        -h|--help)       usage; exit 0 ;;
        *)               echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
    esac
done

if [ -z "${TARGET}" ] && [ -z "${ALL}" ]; then
    echo "ERROR: provide --target NAME or --all." >&2
    usage
    exit 2
fi

cd "${REPO_ROOT}"

if [ -n "${DRY_RUN}" ]; then
    echo "=== DRY RUN — no files will be changed ==="
fi

# Resolve the list of targets to pull.
#
# Kept as a newline-delimited string rather than an array: macOS ships bash
# 3.2, which has no `mapfile` and errors under `set -u` when an empty array is
# expanded.  Target names are PDB-style IDs, so word-splitting is not a risk.
if [ -n "${ALL}" ]; then
    # Ask for directories only, so stray files (README.md, _path_setup.py)
    # are never mistaken for targets.
    TARGET_LIST="$(
        ssh -q -o BatchMode=yes "${HPC_USER_HOST}" \
            "cd ${HPC_BASE}/experiments 2>/dev/null && ls -1d */ 2>/dev/null" \
        | sed 's#/$##' || true
    )"
    if [ -z "${TARGET_LIST}" ]; then
        echo "No targets found under ${HPC_BASE}/experiments on the HPC."
        exit 0
    fi
else
    TARGET_LIST="${TARGET}"
fi

# Read the target list on fd 3, and pass -n to the ssh probe below.  Both are
# needed: commands inside the loop would otherwise consume the target list
# from stdin, and only the first target would ever be synced.
while IFS= read -r t <&3; do
    [ -z "${t}" ] && continue
    REL="experiments/${t}"
    SRC="${HPC_USER_HOST}:${HPC_BASE}/${REL}/"
    DST="${REPO_ROOT}/${REL}/"

    if ! ssh -n -q -o BatchMode=yes "${HPC_USER_HOST}" \
            "[ -d ${HPC_BASE}/${REL} ]" 2>/dev/null; then
        echo "[${t}] SKIP: no remote dir at ${HPC_BASE}/${REL}"
        continue
    fi

    mkdir -p "${DST}"
    echo "[${t}] rsync ${SRC} -> ${DST}"
    # --no-perms/--no-owner/--no-group: the HPC filesystem reports every file
    # as mode 755, so plain `rsync -a` would rewrite the mode of tracked
    # inputs (params.yml, reference PDBs) and show them as modified in git
    # with no content change.  Let the local umask decide instead.
    rsync -av --no-perms --no-owner --no-group ${DRY_RUN} -e ssh \
        --exclude='.DS_Store' \
        --exclude='work/' \
        --exclude='.nextflow*' \
        --exclude='tmp/' \
        "${SRC}" "${DST}"
    echo ""
done 3<<< "${TARGET_LIST}"

echo "Sync from HPC complete."
if [ -n "${DRY_RUN}" ]; then
    echo "(dry-run — nothing was actually pulled)"
fi
