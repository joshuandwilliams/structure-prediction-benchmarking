# containers/

Singularity/Apptainer containers used by the benchmark pipeline.

The built `.img`/`.sif` images are multi-gigabyte and live on the HPC, not
in git (they are gitignored). This directory holds the **definition files**
and build notes needed to recreate them. Add the `.def` files here as they
are written.

## Images the pipeline expects

The image paths are set in [`../nextflow.config`](../nextflow.config) and
can be overridden per-run in the params file:

| Param | Default path (NBI/JIC) | Used by |
|---|---|---|
| `benchmark_container` | `…/Boltz1_Boltz2_Chai1_ColabFold/benchmark_models.img` | Boltz-1, Boltz-2, Chai-1, sequence/constraint/metrics Python |
| `colabfold_container` | `…/ColabFold/colabfold.img` | `COLABFOLD_SEARCH`, ColabFold predictions |

AlphaFold 2-Multimer and AlphaFold 3 are **not** containerised here — they
load from the NBI HPC `source package` system via the `af2_package_id` /
`af3_package_id` params, against the shared reference databases under
`/nbi/Reference-Data/AlphaFold/`.

## Notes

- `benchmark_models.img` bundles Boltz, Chai-1, and the Python stack
  (`gemmi`, `numpy`, `biopython`, `PyYAML`) that the `bin/` scripts need.
  `--no_kernels` is passed to `boltz predict` because `cuequivariance_torch`
  is absent from the image.
- Chai-1 weights are pre-baked at `/opt/chai_cache` inside the image so no
  network download is attempted on the offline compute nodes.

TODO: commit the `.def` files / build recipes for both images.
