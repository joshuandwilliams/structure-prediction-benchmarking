# tests/

Two tiers, mirroring the Mac→HPC development round-trip.

```
tests/
  conftest.py          Shared fixtures: puts bin/ on sys.path; synthetic
                       two-chain PDB builder with known inter-chain distances.
  unit/                Tier 1 — local_unit. Pure-Python unit tests of bin/
                       helpers. Run on every change; no HPC, no GPU.
  characterization/    Tier 2 — hpc. Fixture-based characterization of whole
                       predictor stages against a fresh pipeline run.
  run_pytest.slurm.sh  Submit the suite on the HPC in the lightweight
                       pytest_runner container (pytest + numpy + pandas).
```

## Tier 1 — `local_unit` (Mac, fast)

Pure-Python unit tests of the `bin/` helpers, runnable on a laptop:

- `unit/test_constraint_geometry.py` — Cα parsing, pocket/contact selection,
  YAML block formatting, and both constraint-extractor CLIs end-to-end.
- `unit/test_ipsae.py` — ipSAE / ipAE maths from `compute_metrics.py` (numpy).
- `unit/test_extract_sequences_fasta.py` — the FASTA-mode input adapter.
- `unit/test_validate_boltz_yaml.py` — the pre-flight YAML validator
  (skips cleanly if PyYAML is absent).

```bash
pip install -e '.[test]'      # pytest, numpy, pandas, PyYAML, gemmi, biopython
pytest -m local_unit
```

Some tests `importorskip` an optional parser dependency (`PyYAML`, and —
once added — `gemmi`/`biopython`-backed tests for `extract_sequences.py`),
so the suite stays green on a minimal install and grows coverage when the
full `.[test]` extra is present.

## Tier 2 — `hpc` (cluster, slow)

Characterization tests that pin a predictor stage's parsed metrics against a
committed reference, run against a fresh pipeline run on the HPC (needs the
Singularity containers + GPU). Add fixtures and tests under
`characterization/` as stages are pinned — see
[`characterization/README.md`](characterization/README.md).

```bash
# Runs the suite in the lightweight pytest_runner container (PyYAML tests skip).
sbatch tests/run_pytest.slurm.sh

# A real hpc-tier run that parses predictions with gemmi needs the benchmark
# container instead:
PYTEST_CONTAINER=/hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold/Boltz1_Boltz2_Chai1.img \
    sbatch tests/run_pytest.slurm.sh -m hpc
```

## Adding a unit test

Reuse the `make_pdb_file` / `two_chain_atoms` fixtures in `conftest.py` for
anything that needs coordinates — they give deterministic distances so
expected pocket/contact sets are exact. Keep each `local_unit` test under a
second and free of network/GPU/container dependencies.
