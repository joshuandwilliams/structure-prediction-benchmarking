# containers/

Singularity/Apptainer containers used by the benchmark pipeline.

The built `.img`/`.sif` images are multi-gigabyte and live on the HPC, not
in git (they are gitignored). This directory holds the **definition files**
and build notes needed to recreate them.

- [`Boltz1_Boltz2_Chai1.def`](Boltz1_Boltz2_Chai1.def) — builds the
  `benchmark_container` image (Boltz-1, Boltz-2, Chai-1 + the Python stack).
- [`colabfold.def`](colabfold.def) — builds the `colabfold_container` image
  (MMseqs2 + ColabFold with the AlphaFold2 extra, JAX on the CUDA 12 plugin).
- [`esmfold2.def`](esmfold2.def) — builds the `esmfold2_container` image
  (ESMFold2 on the ESMC-6B backbone, via the `esm` package and the Biohub
  `transformers` fork).

All three are copies of the definition files under
`/hpc-home/jowillia/singularity/`, which is where the images are actually
built. Keep them in step: edit on the HPC, then copy back here.

## Images the pipeline expects

The image paths are set as `params` defaults in [`../main.nf`](../main.nf)
and can be overridden per-run in the params file:

| Param | Default path (NBI/JIC) | Used by |
|---|---|---|
| `benchmark_container` | `…/Boltz1_Boltz2_Chai1_ColabFold/Boltz1_Boltz2_Chai1.img` | Boltz-1, Boltz-2, Chai-1, sequence/constraint/metrics Python |
| `colabfold_container` | `…/ColabFold/colabfold.img` | `COLABFOLD_SEARCH`, ColabFold predictions |
| `esmfold2_container` | `…/ESMFold2/esmfold2.img` | ESMFold2 (`bin/esmfold2_fold.py`); needs a 40 GB+ GPU |

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

Every definition file opens with pre-flight URL checks, so a moved weight URL
aborts the build in seconds rather than after half an hour. All three then
download model weights and pip packages, so they need a build node with
internet access:

```bash
cd /hpc-home/jowillia/singularity/Boltz1_Boltz2_Chai1_ColabFold
singularity build --fakeroot Boltz1_Boltz2_Chai1.img Boltz1_Boltz2_Chai1.def

cd /hpc-home/jowillia/singularity/ColabFold
singularity build --fakeroot colabfold.img colabfold.def

cd /hpc-home/jowillia/singularity/ESMFold2
singularity build --fakeroot esmfold2.img esmfold2.def
```

Approximate build cost: ~13 GB of weights and 30–45 min for the benchmark
image, ~27 GB and longer for ESMFold2 (ESMC-6B alone is ~25 GB).

## Versions in the built images

The definition files pin only some versions explicitly. Boltz, ColabFold and
MMseqs2 are installed unpinned, so the version that actually ran is a property
of the built image rather than of the def file. Read back from the images as
built on 2026-08-15, under Singularity 3.8.7:

| Image | Package | Version |
|---|---|---|
| `Boltz1_Boltz2_Chai1.img` | `boltz` (serves both Boltz-1 and Boltz-2) | 2.2.1 |
| | `chai_lab` | 0.6.1 |
| | `torch` | 2.6.0+cu124 |
| | `trifast` | 0.1.13 |
| | `numpy` / `pandas` / `scipy` | 1.26.4 / 2.3.3 / 1.13.1 |
| | `gemmi` / `biopython` | 0.6.5 / 1.84 |
| `colabfold.img` | `colabfold` | 1.6.1 |
| | `jax` / `jaxlib` | 0.5.3 |
| | MMseqs2 | commit `76da68a` |
| `esmfold2.img` | `esm` | 3.3.0 |
| | `transformers` (Biohub fork) | 4.57.6 |
| | `xformers` | 0.0.29.post3 |
| | `torch` | 2.6.0+cu124 |

AlphaFold2-Multimer and AlphaFold 3 come from the HPC source packages named in
`../main.nf`, at application versions 2.3.2 and 3.0.0 respectively (from each
package's `Singularity.manifest`).

Regenerate this table with:

```bash
singularity exec <image>.img pip list --format=freeze
```
