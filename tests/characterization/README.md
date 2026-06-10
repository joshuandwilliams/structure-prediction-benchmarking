# tests/characterization/  — `hpc` tier

Characterization tests that pin a whole predictor stage's parsed output
against a committed reference, so refactors to `compute_metrics.py` or the
predictor modules can be shown not to change results.

These are **HPC-only** (`@pytest.mark.hpc`): they need real prediction
outputs (and `gemmi`/`numpy` to parse them), which are produced by the
Singularity containers on the GPU queue.

## Intended shape

For each parser in `compute_metrics.py` (`boltz1`, `boltz2`, `chai1`,
`af2m`, `af3`, `colabfold`):

```
characterization/
  <parser>/
    fixtures/                 a small, curated slice of a real prediction
                              directory for that model (a couple of seeds)
    expected_metrics.csv      committed reference output of compute_metrics.py
    test_<parser>_metrics.py  parse fixtures → assert against expected_metrics.csv
```

A test parses the fixture directory with the same `compute_metrics.py`
entrypoint the pipeline uses and asserts the resulting metrics match
`expected_metrics.csv` (within a tolerance for floats). When a change to the
metric code is intentional, regenerate the expected CSV and commit it in the
same change, with the reason in the commit message.

## Status

Not yet populated — the fixtures must be carved from a real benchmark run on
the HPC and pulled back with `scripts/sync_from_hpc.sh`. Until then, the
safety net is the `local_unit` tier (the pure maths and parsing helpers) plus
the byte-for-byte equivalence checks used when refactoring the extractors.
