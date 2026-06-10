# Benchmark target: 6G10_copy

A second run of the **6G10** target (Pik-1 HMA + AVR-PikD, Maqbool et al.
2015), kept in its own experiment directory so its outputs don't collide with
`experiments/6G10/`. Same reference and chain assignment; runs the full
13-model panel including ESMFold2.

| | |
|---|---|
| Input mode | PDB |
| Input | `6G10.pdb` |
| Receptor chain | B |
| Effector chain | C |
| Models | all 13 (see `params.yml`) |

## Provenance

Reference and chain assignment are identical to `experiments/6G10/` — chain B
is the HMA copy bound to AVR-PikD (chain C); chain A is an unbound
crystal-packing copy.

## Outputs

Running `sbatch ../../run_benchmark.slurm.sh ./params.yml` on the HPC writes
`6G10_copy_all_models_results/` here (gitignored). Pull it back for local
analysis with `./scripts/sync_from_hpc.sh --target 6G10_copy` from the repo
root.
