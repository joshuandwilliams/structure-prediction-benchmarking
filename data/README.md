# data/

Reference structures for the benchmark, derived from the NLR structure list in
[`../notes/NLR_Benchmark_PDB_Reference.xlsx`](../notes/NLR_Benchmark_PDB_Reference.xlsx).

```
data/
  all_nlr_pdbs.txt              All 91 PDB IDs from the spreadsheet (tracked).
  benchmark_complexes.tsv       Manifest: the 43 binder–target complexes to
                                extract, with source receptor/target chains
                                (tracked — the source of truth).
  download_solved_structures.sh Pull every PDB in all_nlr_pdbs.txt as mmCIF.
  extract_benchmark_complexes.py Build the 2-chain A/B references from the manifest.

  solved_NLR_structures/        Downloaded mmCIF files (gitignored — bulky,
                                regenerable from the script).
  complexes_for_benchmarking/   The 43 two-chain A/B reference PDBs (tracked).
```

## Workflow

```bash
# 1. Download all 91 structures as mmCIF (needs internet; run on Mac or a login node)
bash data/download_solved_structures.sh

# 2. Build the two-chain benchmark references (needs gemmi: pip install gemmi,
#    or run inside the benchmark container)
python data/extract_benchmark_complexes.py
```

Both scripts are idempotent (skip work already done) and safe to re-run.

## The two chain convention

Every file in `complexes_for_benchmarking/` has **exactly two chains**:

- **chain A = receptor = the PLANT protein** (NLR, integrated HMA domain, or
  host target)
- **chain B = target = the PATHOGEN protein** (effector)

This holds uniformly across all tiers. For Tier-3 effector/host-target
complexes the plant host is still chain A even though biologically the effector
is the "binder" there — uniform `A = plant` labelling keeps every case
trackable and lets a single params template (`receptor_chain: A`,
`effector_chain: B`) drive them all.

`extract_benchmark_complexes.py` chooses, from the chains the manifest lists
for each entity, the receptor/target pair with the most Cα–Cα contacts — so a
multi-copy crystal still yields a genuinely bound pair (e.g. 6G10 picks the
Pikp-1 copy actually bound to AVR-PikD, not the crystal-packing copy).

## The 43 complexes

Tier 1 (18) — small HMA / integrated-domain : effector pairs (Pik, RGA5).
Tier 2 (4)  — other direct NLR-domain : effector binary pairs (Sr35, RPP1, ROQ1, RRS1).
Tier 3 (11) — effector : host-target pairs (plant host = A).
Tier 4 (10) — pairs extracted from larger assemblies (NRC, Sr35/RPP1/WRR4A resistosomes).

See `benchmark_complexes.tsv` for the per-entry tier, system, source chains,
and description. Note `9QT4` and `9QU9` are alternate refinements of the same
WRR4A–CCG40 complex as `9QLV` (kept for completeness; drop them if you want one
representative per complex).

## Updating the list

`all_nlr_pdbs.txt` is generated from the spreadsheet; regenerate it after
editing the sheet. To add a complex to the benchmark set, add a row to
`benchmark_complexes.tsv` (PDB, comma-separated receptor chains, comma-separated
target chains, tier, system, description) and re-run the extractor.
