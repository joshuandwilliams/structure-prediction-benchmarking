# Benchmark target: Pikp1_HMA-AVR_PikF

Pikp-1 HMA x AVR-PikF, run in FASTA mode (no experimental reference complex; confidence metrics only). Receptor = Pikp-1 HMA (entry 1, synthetic chain A), effector = AVR-PikF (entry 2, synthetic chain B).

| | |
|---|---|
| Input mode | FASTA |
| Input | `Pikp1_HMA__AVR_PikF.fasta` |
| Receptor chain | A |
| Effector chain | B |
| Models | see `params.yml` |

## Provenance

TODO: record where the reference came from (PDB accession / how the FASTA
sequences were derived) and any chain-selection rationale.

## Outputs

Running `sbatch ../../run_benchmark.slurm.sh ./params.yml` on the HPC writes
the result tree here (gitignored). Pull it back for local analysis with
`./scripts/sync_from_hpc.sh --target Pikp1_HMA-AVR_PikF` from the repo root.
