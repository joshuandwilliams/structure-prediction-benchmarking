#!/bin/bash
#
# aggregate_seed_outputs.sh
# -------------------------
# Collect a predictor's per-seed output directories into one staging tree.
#
# Each predictor process runs several seeds into sibling directories named
# `output`, `output_seed123`, `output_seed456`, ... This helper, run from the
# process work directory, copies every file with one of the requested
# extensions from each seed directory into `all_outputs/<seed_tag>/<relpath>`,
# preserving the per-seed sub-tree. compute_metrics.py then globs all_outputs/
# recursively across all seeds in a single pass.
#
# This is the shared form of the aggregation loop that used to be duplicated
# verbatim in every Boltz/Chai predictor module.
#
# Usage (from the work dir):
#   aggregate_seed_outputs.sh pdb npz json        # Boltz-1 / Boltz-2
#   aggregate_seed_outputs.sh pdb npz pt json     # Chai-1 (also keeps *.pt)

set -euo pipefail

if [ "$#" -eq 0 ]; then
    echo "Usage: aggregate_seed_outputs.sh EXT [EXT ...]" >&2
    exit 2
fi

# Build a find name-filter expression: \( -name '*.pdb' -o -name '*.npz' ... \)
find_expr=( '(' )
first=1
for ext in "$@"; do
    [ "$first" -eq 0 ] && find_expr+=( -o )
    find_expr+=( -name "*.${ext}" )
    first=0
done
find_expr+=( ')' )

mkdir -p all_outputs
for seed_dir in output output_seed*; do
    [ -d "${seed_dir}" ] || continue
    seed_tag=$(basename "${seed_dir}")
    find "${seed_dir}" "${find_expr[@]}" | while read -r f; do
        rel="${f#${seed_dir}/}"
        dst="all_outputs/${seed_tag}/${rel}"
        mkdir -p "$(dirname "${dst}")"
        cp "${f}" "${dst}"
    done
done
