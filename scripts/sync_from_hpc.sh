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
#   the Mac as authoritative.  Nextflow scratch (work/, .nextflow*) is excluded
#   — only the published *_results/ outputs come back.

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
if [ -n "${ALL}" ]; then
    mapfile -t TARGETS < <(
        ssh -q -o BatchMode=yes "${HPC_USER_HOST}" \
            "ls -1 ${HPC_BASE}/experiments 2>/dev/null" \
        | grep -vE '^(README|_path_setup)' || true
    )
    if [ ${#TARGETS[@]} -eq 0 ]; then
        echo "No targets found under ${HPC_BASE}/experiments on the HPC."
        exit 0
    fi
else
    TARGETS=("${TARGET}")
fi

for t in "${TARGETS[@]}"; do
    REL="experiments/${t}"
    SRC="${HPC_USER_HOST}:${HPC_BASE}/${REL}/"
    DST="${REPO_ROOT}/${REL}/"

    if ! ssh -q -o BatchMode=yes "${HPC_USER_HOST}" \
            "[ -d ${HPC_BASE}/${REL} ]" 2>/dev/null; then
        echo "[${t}] SKIP: no remote dir at ${HPC_BASE}/${REL}"
        continue
    fi

    mkdir -p "${DST}"
    echo "[${t}] rsync ${SRC} -> ${DST}"
    rsync -av ${DRY_RUN} -e ssh \
        --exclude='.DS_Store' \
        --exclude='work/' \
        --exclude='.nextflow*' \
        --exclude='tmp/' \
        "${SRC}" "${DST}"
    echo ""
done

echo "Sync from HPC complete."
if [ -n "${DRY_RUN}" ]; then
    echo "(dry-run — nothing was actually pulled)"
fi
