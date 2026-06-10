# containers/

Singularity/Apptainer containers used by the benchmark pipeline.

The built `.img`/`.sif` images are multi-gigabyte and live on the HPC, not
in git (they are gitignored). This directory holds the **definition files**
and build notes needed to recreate them.

- [`Boltz1_Boltz2_Chai1.def`](Boltz1_Boltz2_Chai1.def) — builds the
  `benchmark_container` image (Boltz-1, Boltz-2, Chai-1 + the Python stack).

## Images the pipeline expects

The image paths are set as `params` defaults in [`../main.nf`](../main.nf)
and can be overridden per-run in the params file:

| Param | Default path (NBI/JIC) | Used by |
|---|---|---|
| `benchmark_container` | `…/Boltz1_Boltz2_Chai1_ColabFold/Boltz1_Boltz2_Chai1.img` | Boltz-1, Boltz-2, Chai-1, sequence/constraint/metrics Python |
| `colabfold_container` | `…/ColabFold/colabfold.img` | `COLABFOLD_SEARCH`, ColabFold predictions |

AlphaFold 2-Multimer and AlphaFold 3 are **not** containerised here — they
load from the NBI HPC `source package` system via the `af2_package_id` /
`af3_package_id` params, against the shared reference databases under
`/nbi/Reference-Data/AlphaFold/`.

## Notes

- `Boltz1_Boltz2_Chai1.img` bundles Boltz, Chai-1, and the Python stack
  (`gemmi`, `numpy`, `biopython`, `PyYAML`) that the `bin/` scripts need.
  `--no_kernels` is passed to `boltz predict` because `cuequivariance_torch`
  is absent from the image.
- Chai-1 weights are pre-baked at `/opt/chai_cache` inside the image so no
  network download is attempted on the offline compute nodes.
- The image name reflects its contents (Boltz-1/Boltz-2/Chai-1); it lives in
  the `…/Boltz1_Boltz2_Chai1_ColabFold/` folder for historical reasons, but
  ColabFold is **not** in it — that's a separate `colabfold.img`. Do not
  confuse it with the sibling `boltz2_negsteer.img` (Boltz-2 only, no Chai-1),
  which belongs to the receptor-resurfacing project.

## Building

The image runs pre-flight URL checks, then downloads ~13 GB of model weights
and pip packages, so it needs internet (build node) and ~30–45 min:

```bash
cd /hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold
singularity build --fakeroot Boltz1_Boltz2_Chai1.img Boltz1_Boltz2_Chai1.def
```

TODO: commit the ColabFold build def (`colabfold.def`) here too.
