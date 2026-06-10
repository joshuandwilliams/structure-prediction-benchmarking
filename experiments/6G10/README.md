# Benchmark target: 6G10

Pik-1 HMA + AVR-PikD complex (Maqbool et al. 2015). Chain B is the HMA copy bound to AVR-PikD (chain C); chain A is an unbound crystal-packing copy.

| | |
|---|---|
| Input mode | PDB |
| Input | `6G10.pdb` |
| Receptor chain | B |
| Effector chain | C |
| Models | see `params.yml` |

## Provenance

TODO: record where the reference came from (PDB accession / how the FASTA
sequences were derived) and any chain-selection rationale.

## Outputs

Running `sbatch ../../run_benchmark.slurm.sh ./params.yml` on the HPC writes
the result tree here (gitignored). Pull it back for local analysis with
`./scripts/sync_from_hpc.sh --target 6G10` from the repo root.
