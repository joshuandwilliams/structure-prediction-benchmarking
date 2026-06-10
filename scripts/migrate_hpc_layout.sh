#!/bin/bash
#
# migrate_hpc_layout.sh  —  ONE-TIME HPC layout migration
# ========================================================
# Brings the existing HPC layout into line with the restructured repo:
#
#   1. Rename the code folder
#         receptor_design/structure-prediction-benchmark
#      -> receptor_design/structure-prediction-benchmarking
#      (preserves the bundled jdk-*/ and nxf_home/ runtime in place).
#
#   2. Move each existing benchmark result tree from
#         receptor_design/benchmarks/<target>/<project>_results/
#      into the repo's gitignored, HPC-side experiment dir
#         receptor_design/structure-prediction-benchmarking/experiments/<target>/
#
# RUN THIS ON THE HPC (it operates on $HOME paths), e.g.:
#       ssh slurm
#       bash ~/receptor_design/structure-prediction-benchmarking/scripts/migrate_hpc_layout.sh --dry
#       bash ~/.../scripts/migrate_hpc_layout.sh        # for real
#
# It is idempotent and guarded: re-running it skips anything already done and
# never overwrites an existing destination. Nextflow scratch (work/,
# .nextflow*) is intentionally LEFT in benchmarks/ — it is regenerable; remove
# it yourself once you're satisfied.
#
# Suggested overall sequence:
#   1. (HPC)  bash scripts/migrate_hpc_layout.sh --dry   # review
#   2. (HPC)  bash scripts/migrate_hpc_layout.sh         # rename + move results
#   3. (Mac)  ./scripts/sync_to_hpc.sh --dry             # review what --delete does
#   4. (Mac)  ./scripts/sync_to_hpc.sh                   # push new code + inputs
#
# Step 3's --delete will NOT touch the migrated results: they live under
# experiments/*/*_results/, which sync_to_hpc.sh excludes.

set -euo pipefail

DRY=0
[ "${1:-}" = "--dry" ] || [ "${1:-}" = "--dry-run" ] && DRY=1
run() { if [ "$DRY" -eq 1 ]; then echo "  [dry] $*"; else echo "  $*"; "$@"; fi; }

BASE="${HOME}/receptor_design"
OLD="${BASE}/structure-prediction-benchmark"
NEW="${BASE}/structure-prediction-benchmarking"
BENCH="${BASE}/benchmarks"

echo "=== 1. Rename code folder ==="
if [ -d "$OLD" ] && [ ! -e "$NEW" ]; then
    run mv "$OLD" "$NEW"
elif [ -d "$NEW" ]; then
    echo "  already at $NEW — skipping rename."
    [ -d "$OLD" ] && echo "  NOTE: $OLD still exists too; investigate before continuing." >&2
else
    echo "ERROR: neither $OLD nor $NEW exists." >&2
    exit 1
fi

echo ""
echo "=== 2. Move benchmark result trees into experiments/<target>/ ==="
if [ ! -d "$BENCH" ]; then
    echo "  no benchmarks dir at $BENCH — nothing to migrate."
    exit 0
fi

shopt -s nullglob
for tdir in "$BENCH"/*/; do
    t="$(basename "$tdir")"
    dest="${NEW}/experiments/${t}"
    echo "[${t}]"
    run mkdir -p "$dest"
    results=( "$tdir"*_results/ )
    if [ ${#results[@]} -eq 0 ]; then
        echo "  (no *_results/ to move)"
        continue
    fi
    for r in "${results[@]}"; do
        rname="$(basename "$r")"
        if [ -e "${dest}/${rname}" ]; then
            echo "  ${rname} already present in experiments/${t} — skipping."
        else
            run mv "$r" "${dest}/${rname}"
        fi
    done
done
shopt -u nullglob

echo ""
if [ "$DRY" -eq 1 ]; then
    echo "Dry run complete — nothing changed."
else
    echo "Migration complete."
    echo "Scratch (work/, .nextflow*, *.out/.err) left under ${BENCH}; remove when ready:"
    echo "    rm -rf ${BENCH}"
fi
