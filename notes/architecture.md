# Architecture notes

A map of how the benchmark pipeline fits together, for orientation before
editing. The authoritative behaviour is always the code in `main.nf` and
`modules/`; this document is the why and the shape.

## Dataflow

```
                    ┌─ EXTRACT_SEQUENCES (PDB mode)
input ──────────────┤                              ──► sequences.json ──► receptor_seq / effector_seq
                    └─ EXTRACT_SEQUENCES_FROM_FASTA (FASTA mode)         (value channels, fanned out)
                                                          │
        ┌─────────────────────────────────────────────────┤
        │ (conditional, only if a selected model needs it) │
   COLABFOLD_SEARCH ──► chain_A.a3m / chain_B.a3m     AF3_SETUP_DB ──► $HOME/af3_db
        │                                                  │
        ▼                                                  ▼
   per-model predictor process  ──►  prediction_dir  ──►  COMPUTE_METRICS (one alias per model)
        (5 seeds × 5 samples)                                   │  ► <model>_metrics.csv
                                                                │  ► <model>_best/
                                                                ▼
                                                         AGGREGATE_RESULTS
                                                          ► all_metrics.csv
                                                          ► all_metrics_ranked_by_effector_rmsd.csv
                                                          ► best_models/
   workflow.onComplete ► predictor_runtime_stats.csv (parsed from trace.txt)
```

## Two input modes

`main.nf` accepts exactly one of `--reference_pdb` (PDB mode) or
`--input_fasta` (FASTA mode). Both preprocessing processes emit an
identically-shaped `sequences.json` **and** a `reference.pdb` — a real
reference in PDB mode, or a REMARK-only placeholder in FASTA mode. Because
`compute_metrics.py` treats a placeholder/empty reference as "no structure
available" and falls through to confidence-only metrics, every downstream
module is mode-agnostic — nothing past preprocessing needs to know which
mode ran.

## Model catalogue & the parser indirection

`ALL_MODELS` lists the 12 runnable model variants; `MODEL_TO_PARSER` maps
each to the **output format** its predictions are in (e.g. `boltz1_msa` and
`boltz1_constrained` both parse as `boltz1`). `COMPUTE_METRICS` is a single
reusable process imported under 12 aliases (`METRICS_BOLTZ1`, …) so each
model publishes to its own `${outdir}/<model>/` subdir while sharing one
`compute_metrics.py` implementation. The `model_tag` (variant name) is
rewritten into the CSV `model` column post-hoc when it differs from the
parser tag, so variants stay distinguishable in the aggregated table.

## Seed loop convention

Every predictor runs 5 seeds (42, 123, 456, 789, 1024) **serially** inside
one process (to respect the small GPU queue), then aggregates all seed
outputs into a single `all_outputs/<seed_tag>/…` staging tree that
`compute_metrics.py` globs recursively. `boltz predict` exits 0 even on
silent parse failures, so each seed is guarded by a "produced ≥1 PDB" check.

The aggregation step is shared: every Boltz/Chai process calls
`bin/aggregate_seed_outputs.sh <ext...>` (Boltz keeps `pdb npz json`; Chai-1
additionally keeps `pt`) instead of carrying its own copy of the loop.

## Constraints (Boltz only)

`extract_constraints_boltz1.py` emits a **pocket-only** block — Boltz-1's
schema forces `max_distance == 6.0` and rejects contact constraints.
`extract_constraints_boltz2.py` emits pocket **+** dense per-residue-pair
contact constraints with arbitrary distances. Both derive restraints from
the reference complex's Cα coordinates; `validate_boltz_yaml.py` re-checks
the generated YAML before `boltz predict` runs.

## Resource accounting

There is no custom stats-collection step. Nextflow's native trace / report /
timeline observers (configured in `nextflow.config`) record per-task
elapsed/cpu/rss/vmem; the `workflow.onComplete` hook in `main.nf` filters
`trace.txt` down to the predictor processes and writes
`predictor_runtime_stats.csv`.

## Known wrinkles

- `main.nf`'s `onComplete` closures (`parse_mem_gb`, `parse_duration`) use
  `return` inside `switch`/loops, which the strict Nextflow language-server
  parser (`nextflow lint`) rejects, though the pipeline runs fine on the
  deployed Nextflow runtime. Treat `nextflow lint` as a regression check
  (don't add *new* errors) rather than a must-be-green gate until those
  closures are rewritten.
- The per-seed `run_seed` shell function is still duplicated across the six
  Boltz processes. Unlike the aggregation loop (now shared via
  `bin/aggregate_seed_outputs.sh`), it wraps the GPU `singularity exec … boltz
  predict` call with per-variant flag differences (`--model boltz1`,
  `--use_potentials`), so extracting it safely needs a real HPC smoke-test to
  verify — deferred until then.
