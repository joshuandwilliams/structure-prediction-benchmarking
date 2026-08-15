# Structure Prediction Benchmarking

A Nextflow DSL2 pipeline that benchmarks protein **complex** structure
predictors against a reference complex. Given a receptor–effector pair
(supplied either as a two-chain reference PDB or as a two-entry FASTA),
it runs each predictor across five seeds, computes confidence and
structural-accuracy metrics, and produces a single ranked table plus the
best model from each predictor.

## Models benchmarked

| Model | Variants | MSA source |
|---|---|---|
| Boltz-1 | baseline, `_msa`, `_constrained` (pocket) | single-seq / shared ColabFold / single-seq |
| Boltz-2 | baseline, `_msa`, `_constrained` (pocket + contacts) | single-seq / shared ColabFold / single-seq |
| Chai-1 | baseline | single-seq |
| AlphaFold2-Multimer | baseline (`af2m`) | full AF2 search |
| AlphaFold 3 | `af3`, `af3_nomsa` | full AF3 search / none |
| ColabFold | `colabfold`, `colabfold_nomsa` | shared ColabFold / none |
| ESMFold2 | `esmfold2` | none (single-sequence) |

The shared ColabFold MSA (`COLABFOLD_SEARCH`) is computed **once** and
reused by every model that needs it (`boltz1_msa`, `boltz2_msa`,
`colabfold`). The AF3 database symlink farm (`AF3_SETUP_DB`) is likewise
built once and shared by `af3` / `af3_nomsa`. GPU jobs are serialized via
`maxForks` to respect the small `jic-gpu` queue at NBI/JIC.

## What the pipeline does

1. **Preprocess** — extract receptor/effector sequences from the
   reference PDB (`EXTRACT_SEQUENCES`) or the input FASTA
   (`EXTRACT_SEQUENCES_FROM_FASTA`) into a shared `sequences.json`.
2. **Shared inputs (conditional)** — run the ColabFold MSA search and/or
   build the AF3 database, only if a selected model needs them.
3. **Predict** — run each selected model (5 seeds × 5 diffusion samples).
4. **Score** — `COMPUTE_METRICS` parses each predictor's output and emits
   confidence metrics (pLDDT, pTM, ipTM, PAE, ipSAE, actifPTM) and, when a
   real reference structure is provided, structural metrics (sequence-aligned
   Cα RMSDs for receptor and effector, including the whole-fit
   effector-on-receptor-aligned RMSD used for ranking).
5. **Aggregate** — concatenate all per-model metrics into `all_metrics.csv`,
   a ranked `all_metrics_ranked_by_effector_rmsd.csv`, and collect each
   model's best structure into `best_models/`. Per-task runtime/memory is
   emitted to `predictor_runtime_stats.csv` from Nextflow's trace.

## Repo layout

```
bin/          Python scripts invoked from Nextflow processes
              (sequence/constraint extraction, YAML validation, metrics)
modules/      Nextflow modules — one file per predictor / stage
experiments/  Per-target benchmark runs. Inputs (params.yml, reference PDB)
              are tracked; run outputs are gitignored and stay on the HPC.
tests/        Unit tests for bin/ helpers, per-model fixture characterization,
              and HPC run wrappers
scripts/      Repo tooling (sync_to_hpc.sh, sync_from_hpc.sh)
containers/   Singularity definition files / build notes (images live on HPC)
notes/        Architecture + design notes
main.nf       Top-level workflow
nextflow.config       Container paths and per-process SLURM resources
params_example.yml        Annotated PDB-mode parameter starter
params_fasta_example.yml  Annotated FASTA-mode parameter starter
run_benchmark.slurm.sh    SLURM submission wrapper
```

## Running the pipeline

The pipeline runs on the HPC; the local Mac repo is for editing only.
Pick an input mode and submit:

```bash
# PDB mode (experimental reference complex → structural + confidence metrics)
sbatch run_benchmark.slurm.sh ./params.yml

# FASTA mode (no reference structure → confidence metrics only)
sbatch run_benchmark.slurm.sh ./params_fasta.yml
```

Copy [`params_example.yml`](params_example.yml) (PDB mode) or
[`params_fasta_example.yml`](params_fasta_example.yml) (FASTA mode) as a
starting point — each documents every parameter inline. Per-target runs
live under [`experiments/`](experiments/README.md).

## The analysis environment

[`environment.yml`](environment.yml) pins the conda environment the Quarto
analyses render in. Versions are exact so a re-render reproduces the
committed figures.

```bash
mamba env create -f environment.yml
mamba activate spb-analysis
cd analysis/structure_prediction_benchmark && quarto render
```

It contains no structure predictor. Those live in the Singularity images
built from [`containers/`](containers/README.md), and AlphaFold2-Multimer and
AlphaFold 3 load from the HPC source-package system.

## Running the tests

Two tiers, mirroring the development round-trip:

- **`pytest -m local_unit`** (Mac, fast). Pure-Python unit tests of the
  `bin/` helpers — ipSAE/RMSD maths, constraint extraction, FASTA/YAML
  parsing. Run on every change during development. Some tests skip cleanly
  if an optional parser dependency (`gemmi`, `PyYAML`) is absent.
- **`pytest -m hpc`** (HPC, slow). Characterization against a fresh
  pipeline run; requires containers + GPU.

```bash
pip install -e '.[test]'   # pytest, numpy, pandas, PyYAML, gemmi, biopython
pytest -m local_unit
```

See [`tests/README.md`](tests/README.md) for the tier layout and how to add
per-model fixture tests.

## Development workflow

Develop on the Mac with Claude Code, sync to the HPC with
`./scripts/sync_to_hpc.sh`, run the pipeline or tests there, iterate. The
HPC has no git; the Mac repo is authoritative and the only place commits
happen. Pull run outputs back for local analysis with
`./scripts/sync_from_hpc.sh --target <NAME>`. Round-trips are slow — batch
related changes and lean on `pytest -m local_unit`, `ruff`, and
`nextflow lint` for fast local feedback.

## License / contact

TODO — license not yet specified.

Maintainer: Joshua Williams (joshuandwilliams@outlook.com).
