#!/bin/bash
#
# Sync the local repo to the HPC via SSH (uses the ~/.ssh/config alias 'slurm').
#
# Usage:
#   ./scripts/sync_to_hpc.sh           # real sync
#   ./scripts/sync_to_hpc.sh --dry     # dry run, no changes
#
# Transport:
#   rsync over SSH, using the 'slurm' host alias defined in ~/.ssh/config
#   (slurm -> slurm.nbi.ac.uk, key id_ed25519_nbi).  Requires passwordless
#   key auth.
#
# Direction & authority:
#   The Mac repo is authoritative.  Develop here, sync up, run the pipeline
#   on the HPC.  The HPC has no git.
#
# Excludes:
#   - .git/, Python/IDE caches, macOS metadata, Office lockfiles
#   - Nextflow runtime bundled on the HPC (jdk-*/, nxf_home/) and per-run
#     scratch (work/, .nextflow*, *_results/) under experiments/<target>/
#   - SLURM logs (*.out, *.err)
#   - Singularity images (*.img, *.sif), multi-GB and HPC-side only
#   - analysis/  The Quarto report is authored and rendered on the Mac from
#     results pulled down by sync_from_hpc.sh.  Nothing on the HPC reads it,
#     and its rendered .html plus the *_files/ JS/CSS bundle are gitignored
#     build artefacts that do not belong on the cluster.  The plotting scripts
#     the HPC does run live in bin/ and are still synced.
#
# The sync uses --delete, so files removed from the local repo are also
# removed on the HPC.  Always run with --dry first when in doubt.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HPC_USER_HOST="slurm"
HPC_DEST="receptor_design/structure-prediction-benchmarking/"

DRY_RUN=""
if [ "${1:-}" = "--dry" ] || [ "${1:-}" = "--dry-run" ]; then
    DRY_RUN="--dry-run"
    echo "=== DRY RUN -- no files will be changed ==="
fi

cd "$REPO_ROOT"

rsync -av --delete $DRY_RUN -e ssh \
    --exclude='.git/' \
    --exclude='__pycache__/' \
    --exclude='.pytest_cache/' \
    --exclude='.ruff_cache/' \
    --exclude='.mypy_cache/' \
    --exclude='*.pyc' \
    --exclude='.DS_Store' \
    --exclude='~$*' \
    --exclude='jdk-*/' \
    --exclude='nxf_home/' \
    --exclude='*.img' \
    --exclude='*.sif' \
    --exclude='data/solved_NLR_structures/' \
    --exclude='analysis/' \
    --exclude='work/' \
    --exclude='.nextflow*' \
    --exclude='*_results/' \
    --exclude='experiments/*/work/' \
    --exclude='experiments/*/tmp/' \
    --exclude='experiments/*/.nextflow*' \
    --exclude='experiments/*/*_results/' \
    --exclude='experiments/*/results/' \
    --exclude='experiments/*/*.out' \
    --exclude='experiments/*/*.err' \
    --exclude='*.out' \
    --exclude='*.err' \
    --exclude='rsync_dryrun*.txt' \
    --exclude='rsync_deletions*.txt' \
    "$REPO_ROOT/" "${HPC_USER_HOST}:${HPC_DEST}"

echo
echo "Sync complete -> ${HPC_USER_HOST}:${HPC_DEST}"
