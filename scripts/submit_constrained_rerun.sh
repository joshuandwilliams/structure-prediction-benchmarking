#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# submit_constrained_rerun.sh
#
# Submit the constrained re-runs (boltz1_constrained + boltz2_constrained)
# prepared by
# setup_constrained_rerun.sh.  Run this ON THE HPC login node.
#
# One SLURM job per target, submitted FROM that target's rerun_b2c/ directory so
# the Nextflow cache and work dir stay inside it (and stay fresh).
#
# Usage:
#   bash scripts/submit_constrained_rerun.sh --dry-run   # print only
#   bash scripts/submit_constrained_rerun.sh             # submit
#   bash scripts/submit_constrained_rerun.sh --only 6G10
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="${REPO_ROOT}/experiments/benchmarks"
SLURM_SCRIPT="${REPO_ROOT}/run_benchmark.slurm.sh"
SUBDIR="rerun_constrained"

DRY=""
ONLY=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run|--dry) DRY="1"; shift ;;
        --only)          ONLY+=("$2"); shift 2 ;;
        -h|--help)       sed -n '2,16p' "$0"; exit 0 ;;
        *)               echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -f "${SLURM_SCRIPT}" ]] || { echo "ERROR: missing ${SLURM_SCRIPT}" >&2; exit 2; }
if [[ -z "${DRY}" ]] && ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch not found — run on the HPC, or use --dry-run." >&2
    exit 2
fi

shopt -s nullglob
n=0
for params in "${BENCH_DIR}"/*/"${SUBDIR}"/params.yml; do
    dir="$(dirname "${params}")"
    pdb="$(basename "$(dirname "${dir}")")"

    if [ ${#ONLY[@]} -gt 0 ]; then
        keep=""
        for o in "${ONLY[@]}"; do [ "$o" = "${pdb}" ] && keep=1; done
        [ -n "${keep}" ] || continue
    fi

    # Refuse to resume onto an existing run: the whole point is a clean cache.
    if [ -d "${dir}/work" ] || [ -d "${dir}/.nextflow" ]; then
        echo "  SKIP ${pdb}: ${SUBDIR}/ already has work/ or .nextflow/." >&2
        echo "       Delete them first if you want a genuinely fresh run." >&2
        continue
    fi

    if [ -n "${DRY}" ]; then
        echo "  (cd ${dir} && sbatch ${SLURM_SCRIPT} params.yml)"
    else
        jid="$(cd "${dir}" && sbatch --parsable "${SLURM_SCRIPT}" params.yml)"
        echo "  ${pdb}: job ${jid}"
    fi
    n=$((n + 1))
done

echo
if [ -n "${DRY}" ]; then
    echo "DRY RUN — ${n} job(s) would be submitted, nothing was."
else
    echo "Submitted ${n} job(s)."
    cat <<'EOF'

Watch:
  squeue -u $USER
  tail -f experiments/benchmarks/<PDB>/rerun_constrained/.nextflow.log

When they finish, check whether the constraints actually did anything:
  bash scripts/verify_constrained_rerun.sh
EOF
fi
