# Benchmark target: 7QZD

Pik-1 HMA / AVR-Pik complex. Receptor (HMA) = chain A, effector = chain C.

| | |
|---|---|
| Input mode | PDB |
| Input | `7QZD.pdb` |
| Receptor chain | A |
| Effector chain | C |
| Models | see `params.yml` |

## Provenance

TODO: record where the reference came from (PDB accession / how the FASTA
sequences were derived) and any chain-selection rationale.

## Outputs

Running `sbatch ../../run_benchmark.slurm.sh ./params.yml` on the HPC writes
the result tree here (gitignored). Pull it back for local analysis with
`./scripts/sync_from_hpc.sh --target 7QZD` from the repo root.
