#!/bin/bash
# Regenerate every committed CSV under data/metrics/.
#
# These are what the Quarto analyses read, committed so they render without a
# pipeline run. All but the last two are derived from benchmark run outputs, so
# they are only valid for the run that produced them and must be rebuilt after
# any re-run.
#
# Run from anywhere. Needs the spb-analysis environment, plus EMBOSS on PATH
# for the reference-similarity step and network access for release dates.
#
#   bash scripts/build_metrics_csvs.sh              # all of them
#   bash scripts/build_metrics_csvs.sh --skip-runs  # only the run-independent two

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${REPO}/bin"
OUT="${REPO}/data/metrics"
BENCH="${REPO}/experiments/benchmarks"

SKIP_RUNS=0
[ "${1:-}" = "--skip-runs" ] && SKIP_RUNS=1

mkdir -p "${OUT}"

if [ "${SKIP_RUNS}" -eq 0 ]; then
    if [ ! -d "${BENCH}" ]; then
        echo "ERROR: ${BENCH} not found. Run the benchmark first, or pass --skip-runs." >&2
        exit 1
    fi

    # One row per (model, msa, pdb): the most confident of its 25 predictions.
    python "${BIN}/combine_metrics.py" \
        --benchmarks-dir "${BENCH}" \
        --output "${OUT}/combined_metrics.csv"

    # Every prediction, unselected, for threshold sweeps and cost per prediction.
    python "${BIN}/build_per_prediction_csv.py" \
        --benchmarks-dir "${BENCH}" \
        --output "${OUT}/per_prediction_metrics.csv"

    # Wall-clock and peak memory, parsed from each run's Nextflow trace.
    python "${BIN}/build_runtime_csv.py" \
        --benchmarks-dir "${BENCH}" \
        --output "${OUT}/runtime_stats.csv"

    # Which receptor surface each confidence-selected pose landed on.
    # The output name differs from the script default, which is why this
    # script exists rather than the commands living in shell history.
    python "${BIN}/classify_predicted_interface.py" \
        --output "${OUT}/predicted_interfaces_plddt.csv"

    # Whether Boltz-2's no-MSA failures collapse onto the AVR-Pia surface.
    python "${BIN}/common_wrong_interface.py" \
        --model boltz2 --msa no_msa --frame 6Q76 --system Pik \
        --output "${OUT}/common_interface_boltz2_no_msa.csv"
fi

# Independent of any run: built from the reference PDBs and the RCSB API.
python "${BIN}/compute_reference_similarity.py" \
    --output "${OUT}/reference_similarity.csv"

# Run only to refresh the cached release dates. Its plot is a side effect and
# is discarded, so data/metrics/ stays data.
python "${BIN}/plot_training_set_membership.py" \
    --dates-csv "${OUT}/pdb_release_dates.csv" \
    --outdir "$(mktemp -d)"

echo
echo "Wrote:"
ls -1 "${OUT}"
