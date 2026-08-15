#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# verify_constrained_rerun.sh
#
# Answer the question the re-run exists to answer: with correctly numbered
# constraints, do they change anything?
#
# For each target and each constrained variant, compare the re-run's published
# structure against the ORIGINAL UNCONSTRAINED structure of the same base model
# from the main benchmark.
#
#   DIFFERENT  -> the run genuinely re-executed.  Necessary, not sufficient:
#                 Boltz is stochastic, so a fresh unconstrained run would differ
#                 too.  To show the constraints HELPED, compare ra_eff in
#                 combined_metrics.csv once the results are pulled back.
#   IDENTICAL  -> decisive failure.  Byte-identical output from a different
#                 input cannot happen by chance, so the constraints are still
#                 being discarded and the next step is the Boltz constraint
#                 schema, not another re-run.
#
# Run wherever both result trees exist (HPC, or locally after a pull).
#
# Usage:  bash scripts/verify_constrained_rerun.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${REPO_ROOT}/experiments/benchmarks"
SUBDIR="rerun_constrained"

md5of() {
    if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'
    else md5 -q "$1"; fi
}

printf "%-8s %-10s %-12s %s\n" "TARGET" "MODEL" "VERDICT" "DETAIL"
same=0; diff=0; missing=0

shopt -s nullglob
for params in "${BENCH_DIR}"/*/"${SUBDIR}"/params.yml; do
    dir="$(dirname "${params}")"
    pdb="$(basename "$(dirname "${dir}")")"

    for base in boltz1 boltz2; do
        new="$(ls "${dir}/${pdb}_constrained_rerun_results/best_models/${base}_constrained_best/"*.pdb 2>/dev/null | head -1 || true)"
        old="$(ls "${BENCH_DIR}/${pdb}/${pdb}_benchmark_results/best_models/${base}_best/"*.pdb 2>/dev/null | head -1 || true)"

        if [ -z "${new}" ]; then
            printf "%-8s %-10s %-12s %s\n" "${pdb}" "${base}" "NO-RESULT" "re-run produced nothing yet"
            missing=$((missing + 1)); continue
        fi
        if [ -z "${old}" ]; then
            printf "%-8s %-10s %-12s %s\n" "${pdb}" "${base}" "NO-BASELINE" "original unconstrained not found"
            missing=$((missing + 1)); continue
        fi

        if [ "$(md5of "${new}")" = "$(md5of "${old}")" ]; then
            printf "%-8s %-10s %-12s %s\n" "${pdb}" "${base}" "IDENTICAL" "constraints STILL discarded"
            same=$((same + 1))
        else
            printf "%-8s %-10s %-12s %s\n" "${pdb}" "${base}" "DIFFERENT" "run re-executed"
            diff=$((diff + 1))
        fi
    done
done

cat <<EOF

  DIFFERENT: ${diff}    IDENTICAL: ${same}    missing: ${missing}

Any IDENTICAL row means that variant is still discarding its constraints.

If all rows are DIFFERENT, that only confirms the runs executed. The question
of whether the constraints IMPROVED anything is answered by comparing ra_eff
between the constrained re-run and the unconstrained baseline, which needs the
results pulling back:

    ./scripts/sync_from_hpc.sh --all
EOF
