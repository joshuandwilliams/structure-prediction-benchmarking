#!/bin/bash
#SBATCH --job-name="rmsd_6G10"
#SBATCH -p jic-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=rmsd_6G10_%j.out
#SBATCH --error=rmsd_6G10_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jowillia@nbi.ac.uk

# =============================================================================
# Re-run compute_metrics.py against every model's predictions, this time
# WITH --reference-pdb 6G10.pdb, then concatenate all per-model CSVs and
# rank by rmsd_effector_receptor_aligned (ascending — lowest first).
#
# Why this works without re-running the pipeline:
#   Each METRICS_<MODEL> Nextflow process stages its upstream prediction
#   directory in as ./predictions/, so every PDB/CIF the original metrics
#   step saw is still reachable from those work dirs (via symlink). The
#   pdb_path column in the original CSV is literally relative to that
#   work dir. We just call compute_metrics.py again from inside each
#   metrics work dir with --reference-pdb pointing at 6G10.pdb.
#
# Place 6G10.pdb in the directory you sbatch this script from.
# =============================================================================

set -euo pipefail

# ── Paths ────────────────────────────────────────────────────────────────────
PIPELINE_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking"
# NOTE: this was a one-off analysis of the original Pikp1 benchmark run; its
# outputs live alongside this script in rmsd_vs_6G10_results/. To re-run, point
# WORK_ROOT at a run's work/ tree and refresh the per-model work-dir hashes in
# the MODELS array below from that run's .nextflow.log (grep recipe below).
BENCHMARK_DIR="/hpc-home/jowillia/receptor_design/structure-prediction-benchmarking/experiments/Pikp1_HMA-AVR_PikF"
WORK_ROOT="${BENCHMARK_DIR}/work"

# Read benchmark_container from main.nf, where it's set as a params default.
# The params YAML doesn't override it for this experiment, so main.nf is the
# source of truth — and reading it dynamically means this script tracks any
# future container path change in the pipeline automatically.
BENCHMARK_CONTAINER="$(awk '
    /^[[:space:]]*params\.benchmark_container[[:space:]]*=/ {
        sub(/^[^=]*=[[:space:]]*/, "")
        gsub(/^["'\'']|["'\'']$/, "")
        sub(/[[:space:]]*\/\/.*$/, "")
        gsub(/["'\'']/, "")
        print
        exit
    }
' "${PIPELINE_DIR}/main.nf")"

if [ -z "${BENCHMARK_CONTAINER}" ]; then
    echo "ERROR: could not extract params.benchmark_container from ${PIPELINE_DIR}/main.nf"
    exit 1
fi
if [ ! -f "${BENCHMARK_CONTAINER}" ]; then
    echo "ERROR: benchmark container not found at ${BENCHMARK_CONTAINER}"
    exit 1
fi

REFERENCE_PDB="$(realpath ./6G10.pdb)"
REF_RECEPTOR_CHAIN="B"   # Pikp1-HMA   (chain B in 6G10; chain A is unbound copy)
REF_EFFECTOR_CHAIN="C"   # AVR-PikF    (chain C in 6G10)

# The benchmark predictions use synthetic chain labels A (receptor) and B
# (effector) from FASTA-mode. compute_metrics.py uses --receptor-chain /
# --effector-chain to look up BOTH the predicted and reference chains under
# the same ID, so we'll briefly remap 6G10's B→A and C→B into a temp PDB.
PRED_RECEPTOR_CHAIN="A"
PRED_EFFECTOR_CHAIN="B"

# Output dir — must be absolute because we pushd into per-model workdirs
OUT_DIR="$(realpath -m ./rmsd_vs_6G10_results)"
PER_MODEL_DIR="${OUT_DIR}/per_model"
mkdir -p "${PER_MODEL_DIR}"

# ── Sanity checks (need 6G10 present before parsing it) ──────────────────────
if [ ! -f "${REFERENCE_PDB}" ]; then
    echo "ERROR: 6G10.pdb not found at ${REFERENCE_PDB}"
    echo "Place 6G10.pdb in the directory you sbatch this script from."
    exit 1
fi
if [ ! -f "${PIPELINE_DIR}/bin/compute_metrics.py" ]; then
    echo "ERROR: compute_metrics.py not found at ${PIPELINE_DIR}/bin/"
    exit 1
fi

# ── Extract chain lengths (Cα count) from 6G10 ───────────────────────────────
# These are only used for ipSAE/ipae bookkeeping, not structural RMSD, but
# we read them from the reference PDB so they're always correct rather than
# hardcoded.  Counts unique residues (chain + resSeq + iCode) with a CA atom
# in the named chain, ignoring HETATM.
get_chain_len() {
    local pdb="$1" chain="$2"
    awk -v ch="$chain" '
        /^ATOM/ && substr($0, 13, 4) == " CA " && substr($0, 22, 1) == ch {
            key = substr($0, 22, 1) "_" substr($0, 23, 5)
            if (!(key in seen)) { seen[key] = 1; n++ }
        }
        END { print n+0 }
    ' "$pdb"
}

CHAIN_A_LEN="$(get_chain_len "${REFERENCE_PDB}" "${REF_RECEPTOR_CHAIN}")"
CHAIN_B_LEN="$(get_chain_len "${REFERENCE_PDB}" "${REF_EFFECTOR_CHAIN}")"

if [ "${CHAIN_A_LEN}" -eq 0 ] || [ "${CHAIN_B_LEN}" -eq 0 ]; then
    echo "ERROR: failed to extract chain lengths from ${REFERENCE_PDB}"
    echo "  receptor chain ${REF_RECEPTOR_CHAIN}: ${CHAIN_A_LEN} CA atoms"
    echo "  effector chain ${REF_EFFECTOR_CHAIN}: ${CHAIN_B_LEN} CA atoms"
    exit 1
fi

# ── Build a remapped reference PDB ───────────────────────────────────────────
# compute_metrics.py uses --receptor-chain / --effector-chain to look up the
# SAME chain IDs in both predicted and reference structures. The predictions
# use A/B but 6G10 uses B/C (with an unbound A copy of the receptor we want
# to discard). Write a temp PDB containing only 6G10 chain B → A and C → B.
REMAPPED_REF="${OUT_DIR}/6G10_remapped_AB.pdb"
mkdir -p "$(dirname "${REMAPPED_REF}")"
awk -v rec="${REF_RECEPTOR_CHAIN}" -v eff="${REF_EFFECTOR_CHAIN}" '
    /^ATOM/ {
        c = substr($0, 22, 1)
        if (c == rec) { print substr($0,1,21) "A" substr($0,23) }
        else if (c == eff) { print substr($0,1,21) "B" substr($0,23) }
        next
    }
    /^TER/ {
        c = substr($0, 22, 1)
        if (c == rec) { print substr($0,1,21) "A" substr($0,23) }
        else if (c == eff) { print substr($0,1,21) "B" substr($0,23) }
        next
    }
    /^END/ { print }
' "${REFERENCE_PDB}" > "${REMAPPED_REF}"

# Sanity check: confirm the remapped file has Cα counts in A and B that
# match what we extracted from B and C of the original.
REMAP_A_LEN="$(get_chain_len "${REMAPPED_REF}" "A")"
REMAP_B_LEN="$(get_chain_len "${REMAPPED_REF}" "B")"
if [ "${REMAP_A_LEN}" -ne "${CHAIN_A_LEN}" ] || [ "${REMAP_B_LEN}" -ne "${CHAIN_B_LEN}" ]; then
    echo "ERROR: remapped reference chain count mismatch"
    echo "  expected A=${CHAIN_A_LEN} B=${CHAIN_B_LEN}"
    echo "  got      A=${REMAP_A_LEN} B=${REMAP_B_LEN}"
    exit 1
fi

# ── METRICS_<MODEL> work directories  (model_tag → parser_tag → workdir) ────
# Pulled from .nextflow.log for the run that produced
# all_metrics_ranked_by_effector_rmsd.csv. If you re-run the pipeline, refresh
# these by grepping the new log:
#
#   grep -oE 'METRICS_[A-Z0-9_]+ \([a-z0-9_]+\).*workDir: [^ ]+' .nextflow.log \
#       | sort -u
#
# Format: "model_tag|parser_tag|workdir"
MODELS=(
  "af2m|af2m|${WORK_ROOT}/97/408a63b04be0f4753427f1d42691f1"
  "af3|af3|${WORK_ROOT}/c3/ea8301095334895f55e89d669754a9"
  "af3_nomsa|af3|${WORK_ROOT}/cc/4783ee6802514db4cda95060eb3f07"
  "boltz1|boltz1|${WORK_ROOT}/bd/1c5f785e129d81b95c592d8bccb3a9"
  "boltz1_msa|boltz1|${WORK_ROOT}/b5/9a09d79e681506e05939d9d5e97f48"
  "boltz2|boltz2|${WORK_ROOT}/35/2ba127983178d3f796f1ff2bc83ed3"
  "boltz2_msa|boltz2|${WORK_ROOT}/f5/1f593d6067dd6b10a10b9966da2e0d"
  "chai1|chai1|${WORK_ROOT}/7f/20b0522be03504f9d614f28c528ac0"
  "colabfold|colabfold|${WORK_ROOT}/dc/1113a032847bcef110a703b943ebcf"
  "colabfold_nomsa|colabfold|${WORK_ROOT}/cf/2cd1a36af05b42359c6d270a07f5d5"
)

echo "============================================================"
echo "RMSD vs 6G10 — receptor-aligned effector ranking"
echo "============================================================"
echo "Reference PDB:    ${REFERENCE_PDB}"
echo "  receptor:       chain ${REF_RECEPTOR_CHAIN} (${CHAIN_A_LEN} residues)"
echo "  effector:       chain ${REF_EFFECTOR_CHAIN} (${CHAIN_B_LEN} residues)"
echo "Remapped ref:     ${REMAPPED_REF}"
echo "  receptor → A, effector → B  (matches predictions)"
echo "Container:        ${BENCHMARK_CONTAINER}"
echo "Output dir:       ${OUT_DIR}"
echo "Date:             $(date)"
echo "============================================================"
echo ""

# ── Run compute_metrics.py against each model's predictions/ tree ───────────
for entry in "${MODELS[@]}"; do
    IFS='|' read -r MODEL_TAG PARSER_TAG WORKDIR <<< "${entry}"

    echo "── ${MODEL_TAG} ────────────────────────────────────────────"
    echo "  workdir:    ${WORKDIR}"

    if [ ! -d "${WORKDIR}" ]; then
        echo "  SKIP: work dir does not exist"
        echo ""
        continue
    fi
    if [ ! -d "${WORKDIR}/predictions" ]; then
        echo "  SKIP: ${WORKDIR}/predictions/ not found"
        echo ""
        continue
    fi

    OUT_CSV="${PER_MODEL_DIR}/${MODEL_TAG}_metrics_vs_6G10.csv"

    # Run compute_metrics.py *from inside the metrics work dir* so that the
    # pdb_path values it writes match the original CSV's relative paths
    # (relative to the metrics workdir, where predictions/ is the staged tree).
    pushd "${WORKDIR}" > /dev/null

    # Bind both the work dir (for predictions/) and the launch dir (for
    # the reference PDB and the output CSV).
    singularity exec \
        --bind "${WORKDIR}:${WORKDIR}" \
        --bind "$(dirname "${REMAPPED_REF}"):$(dirname "${REMAPPED_REF}")" \
        --bind "${PER_MODEL_DIR}:${PER_MODEL_DIR}" \
        "${BENCHMARK_CONTAINER}" \
        python "${PIPELINE_DIR}/bin/compute_metrics.py" \
            --model "${PARSER_TAG}" \
            --prediction-dir predictions \
            --chain-lengths "${CHAIN_A_LEN}" "${CHAIN_B_LEN}" \
            --output-csv "${OUT_CSV}" \
            --reference-pdb "${REMAPPED_REF}" \
            --receptor-chain "${PRED_RECEPTOR_CHAIN}" \
            --effector-chain "${PRED_EFFECTOR_CHAIN}" \
        || { echo "  ERROR: compute_metrics.py failed for ${MODEL_TAG}"; popd > /dev/null; echo ""; continue; }

    popd > /dev/null

    # Rewrite the 'model' column for variant tags (boltz1_msa, af3_nomsa, etc.)
    # Same trick metrics.nf uses.
    if [ "${PARSER_TAG}" != "${MODEL_TAG}" ]; then
        python3 - "${OUT_CSV}" "${MODEL_TAG}" << 'PYEOF'
import csv, sys
csv_path, tag = sys.argv[1], sys.argv[2]
with open(csv_path) as f:
    rows = list(csv.reader(f))
if not rows: sys.exit(0)
header = rows[0]
mcol = header.index("model") if "model" in header else 0
for row in rows[1:]:
    if row: row[mcol] = tag
with open(csv_path, "w", newline="") as f:
    csv.writer(f).writerows(rows)
PYEOF
    fi

    # Annotate each row with the metrics work dir it came from, so we can
    # resolve symlinks to absolute paths in the merge step.
    python3 - "${OUT_CSV}" "${WORKDIR}" << 'PYEOF'
import csv, sys
csv_path, workdir = sys.argv[1], sys.argv[2]
with open(csv_path) as f:
    rows = list(csv.reader(f))
if not rows: sys.exit(0)
header = rows[0]
if "metrics_workdir" not in header:
    header.append("metrics_workdir")
    for row in rows[1:]:
        row.append(workdir)
with open(csv_path, "w", newline="") as f:
    csv.writer(f).writerows(rows)
PYEOF

    n_rows=$(($(wc -l < "${OUT_CSV}") - 1))
    echo "  wrote ${n_rows} rows → ${OUT_CSV}"
    echo ""
done

# ── Merge, resolve absolute paths, rank ─────────────────────────────────────
echo "============================================================"
echo "Merging per-model CSVs and resolving absolute PDB paths..."
echo "============================================================"

MERGED_CSV="${OUT_DIR}/all_metrics_vs_6G10.csv"
RANKED_CSV="${OUT_DIR}/all_metrics_ranked_by_effector_rmsd_vs_6G10.csv"
PATH_LIST="${OUT_DIR}/ranked_pdb_paths.txt"

python3 - "${PER_MODEL_DIR}" "${MERGED_CSV}" "${RANKED_CSV}" "${PATH_LIST}" << 'PYEOF'
import csv, glob, os, sys

per_model_dir, merged_csv, ranked_csv, path_list = sys.argv[1:5]

csv_files = sorted(glob.glob(os.path.join(per_model_dir, "*_metrics_vs_6G10.csv")))
if not csv_files:
    print("ERROR: no per-model CSVs found")
    sys.exit(1)

all_rows = []
header = None
for cf in csv_files:
    with open(cf) as f:
        reader = csv.reader(f)
        h = next(reader, None)
        if h is None:
            continue
        if header is None:
            header = h
        for row in reader:
            if not row:
                continue
            # Pad short rows to header length
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))
            all_rows.append(row)

if header is None or not all_rows:
    print("ERROR: no rows collected")
    sys.exit(1)

# Resolve absolute path: pdb_path is relative to metrics_workdir.
# Follow symlinks all the way to the real file under work/<hash>/...
pdb_idx     = header.index("pdb_path")
wdir_idx    = header.index("metrics_workdir")
rmsd_idx    = header.index("rmsd_effector_receptor_aligned")
model_idx   = header.index("model")
mname_idx   = header.index("model_name")

if "abs_pdb_path" not in header:
    header.append("abs_pdb_path")
    for row in all_rows:
        rel = row[pdb_idx]
        wdir = row[wdir_idx]
        cand = os.path.join(wdir, rel) if rel else ""
        if cand and os.path.exists(cand):
            row.append(os.path.realpath(cand))
        else:
            row.append("MISSING")

# Merged (unranked, all rows kept including any with empty rmsd)
with open(merged_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(all_rows)

# Rank by rmsd_effector_receptor_aligned ASC. Drop rows where it's empty/None
# (those are rows compute_metrics.py couldn't compute an RMSD for, e.g. parser
# fallbacks). Sort numerically.
def rmsd_key(row):
    v = row[rmsd_idx]
    try:
        return float(v)
    except (ValueError, TypeError):
        return float("inf")

ranked = sorted(all_rows, key=rmsd_key)
ranked = [r for r in ranked if rmsd_key(r) != float("inf")]

with open(ranked_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    w.writerows(ranked)

# Plain text list: rank, model, model_name, rmsd, abs_path
abs_idx = header.index("abs_pdb_path")
with open(path_list, "w") as f:
    f.write(f"{'rank':>5}  {'model':<18} {'model_name':<32} {'eff_rmsd':>10}  abs_pdb_path\n")
    f.write("-" * 120 + "\n")
    for i, row in enumerate(ranked, 1):
        f.write(f"{i:>5}  {row[model_idx]:<18} {row[mname_idx]:<32} "
                f"{row[rmsd_idx]:>10}  {row[abs_idx]}\n")

print(f"Merged:  {merged_csv}  ({len(all_rows)} rows)")
print(f"Ranked:  {ranked_csv}  ({len(ranked)} rows with valid RMSD)")
print(f"Paths:   {path_list}")
print()
print("Top 20 by receptor-aligned effector RMSD vs 6G10:")
print(f"{'rank':>5}  {'model':<18} {'model_name':<32} {'eff_rmsd':>10}")
print("-" * 75)
for i, row in enumerate(ranked[:20], 1):
    print(f"{i:>5}  {row[model_idx]:<18} {row[mname_idx]:<32} {row[rmsd_idx]:>10}")
PYEOF

echo ""
echo "============================================================"
echo "Done: $(date)"
echo "============================================================"
echo "Outputs:"
echo "  ${OUT_DIR}/all_metrics_vs_6G10.csv"
echo "  ${OUT_DIR}/all_metrics_ranked_by_effector_rmsd_vs_6G10.csv"
echo "  ${OUT_DIR}/ranked_pdb_paths.txt"
echo "  ${OUT_DIR}/per_model/*.csv"
