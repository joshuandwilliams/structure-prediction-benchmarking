#!/bin/bash
#
# download_solved_structures.sh
# -----------------------------
# Download every structure in the NLR benchmark list (data/all_nlr_pdbs.txt)
# into data/solved_NLR_structures/ as mmCIF (.cif).
#
# Why mmCIF and not legacy .pdb:
#   - It is universal: every PDB entry has an mmCIF file, including the large
#     cryo-EM assemblies (resistosomes) that have NO legacy .pdb at all.
#   - Its chain naming (auth_asym_id) matches data/benchmark_complexes.tsv,
#     so extract_benchmark_complexes.py can find the requested chains. (Legacy
#     .pdb files sometimes rename chains, e.g. 7A8W's AAA/CCC.)
# The 2-chain benchmark references are written back out as .pdb by the
# extractor — that is the format the pipeline's constraint parser reads.
#
# Skips anything already present, so it is safe to re-run. Needs internet —
# run on the Mac (or an HPC login node), then sync the data/ folder up.

set -euo pipefail

DATA="$(cd "$(dirname "${BASH_SOURCE[0]}")/../data" && pwd)"
LIST="${DATA}/all_nlr_pdbs.txt"
OUT="${DATA}/solved_NLR_structures"
BASE="https://files.rcsb.org/download"

mkdir -p "${OUT}"

n=0; ok=0; fail=0
while IFS= read -r raw; do
    id="$(printf '%s' "${raw}" | tr -d '[:space:]' | tr 'a-z' 'A-Z')"
    [ -z "${id}" ] && continue
    case "${id}" in \#*) continue ;; esac
    n=$((n + 1))

    if [ -s "${OUT}/${id}.cif" ]; then
        echo "[${id}] already present"; ok=$((ok + 1)); continue
    fi

    if curl -fsSL "${BASE}/${id}.cif" -o "${OUT}/${id}.cif" 2>/dev/null && [ -s "${OUT}/${id}.cif" ]; then
        echo "[${id}] cif"; ok=$((ok + 1))
    else
        rm -f "${OUT}/${id}.cif"
        echo "[${id}] FAILED to download" >&2; fail=$((fail + 1))
    fi
done < "${LIST}"

echo ""
echo "Done: ${ok}/${n} present, ${fail} failed -> ${OUT}"
