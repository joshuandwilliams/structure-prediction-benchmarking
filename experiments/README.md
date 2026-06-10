# experiments/

Working space for the benchmark runs themselves. The production pipeline
(top-level `main.nf` + `modules/`) is a fixed transformation from a
receptor–effector input to a ranked metrics table; this directory is where
that transformation is *applied* to specific targets, and where the outputs
are kept for analysis.

One subdirectory per **target** — a reference complex (or FASTA pair) being
benchmarked across the model panel.

## Directory layout

```
experiments/
  _path_setup.py     Makes the repo's bin/ importable from analysis scripts
                     (so they can reuse compute_metrics helpers instead of
                     reimplementing ipSAE / RMSD maths).

  <target>/          One directory per benchmarked complex.
    README.md        What this target is, chain assignment, provenance.
    params.yml       The parameter file passed to the pipeline for this run.
                     Paths are relative to this directory.
    <target>.pdb     Reference complex (PDB mode), tracked — small.
    <target>.fasta   Input FASTA (FASTA mode), if applicable.

    <project>_results/   Pipeline output tree. GITIGNORED — produced on the
                         HPC, pulled back with scripts/sync_from_hpc.sh for
                         local analysis, never committed (large).
    work/  .nextflow*    Nextflow scratch. GITIGNORED, HPC-side only.
```

## What is tracked vs. gitignored

**Tracked** (small, authoritative, lives on the Mac and syncs up):
the per-target `params.yml`, the reference PDB / input FASTA, and the
`README.md`.

**Gitignored** (large, produced on the HPC, pulled back only for analysis):
everything matching `experiments/*/*_results/`, `experiments/*/work/`,
`experiments/*/.nextflow*`, and SLURM logs. See the repo `.gitignore`.

## Running a target

From the HPC, inside the target directory:

```bash
cd experiments/<target>
sbatch ../../run_benchmark.slurm.sh ./params.yml
```

The run writes `<project>_results/` here. Pull it back to the Mac for
analysis with:

```bash
./scripts/sync_from_hpc.sh --target <target>
```

## Mac and HPC

The Mac copy is authoritative for the tracked inputs; the HPC produces the
results. `scripts/sync_to_hpc.sh` pushes inputs up (excluding result trees);
`scripts/sync_from_hpc.sh` pulls a target's results down. Run outputs are
never committed.

## Analysis scripts

Analysis/plotting scripts that compare predictors across a target's results
should reuse the metric logic in `bin/` rather than reimplementing it.
Import the path shim once at the top of a script:

```python
import _path_setup  # noqa: F401  — adds repo bin/ to sys.path
from compute_metrics import compute_ipsae   # etc.
```
