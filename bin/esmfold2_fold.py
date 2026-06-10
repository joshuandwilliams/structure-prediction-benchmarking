#!/usr/bin/env python3
"""
esmfold2_fold.py
----------------
Run one seeded ESMFold2 fold of a two-chain (receptor + effector) complex.

ESMFold2 (CZ Biohub, ESMC-6B backbone) is a single-sequence, MSA-free,
diffusion-based *complex* predictor. There is no MSA input — the two chains
are passed as two ProteinInput entries and folded directly. This script is
invoked inside esmfold2.img (which provides the `esm` package and the Biohub
`transformers` fork) once per seed by modules/esmfold2.nf.

It reads the receptor/effector sequences from a small JSON file and writes,
into --out-dir:

    esmfold2_pred.cif   the predicted complex (mmCIF)
    confidences.json    {"plddt": [...0-100], "ptm": float, "iptm"?, "pae"?}

compute_metrics.py's parse_esmfold2 consumes that pair, exactly as it does the
AF3 model.cif + confidences.json.

Usage:
    python esmfold2_fold.py --input-json input.json --out-dir output_seed42 \
        --seed 42 [--num-loops 20] [--num-sampling-steps 100] [--repo biohub/ESMFold2]

input.json shape:
    {"receptor": {"id": "A", "sequence": "MKFL..."},
     "effector": {"id": "B", "sequence": "GTAL..."}}
"""

from __future__ import annotations

import argparse
import json
import os
import sys


def _to_float_list(x):
    """Coerce a torch tensor / numpy array / scalar / list to a flat float list."""
    import numpy as np

    if hasattr(x, "detach"):          # torch tensor
        x = x.detach().cpu().numpy()
    return np.asarray(x, dtype=float).ravel().tolist()


def build_confidences(plddt, ptm=None, iptm=None, pae=None):
    """Assemble the confidences dict written next to the structure.

    pLDDT is rescaled from ESMFold2's 0-1 range to the 0-100 scale the rest of
    the pipeline (and compute_metrics.py) expects; pTM/ipTM are kept on their
    native 0-1 scale. ``pae``, when the model exposes it, is emitted as a
    nested list so ipSAE/ipAE can be computed downstream — when it is absent
    the interface metrics simply fall back to 0, as for any model without a
    PAE matrix.

    This is the pure, GPU-free part of the runner and is unit-tested directly.
    """
    import numpy as np

    vals = _to_float_list(plddt)
    if vals and max(vals) <= 1.5:          # looks like a 0-1 fraction → 0-100
        vals = [v * 100.0 for v in vals]
    conf = {"plddt": vals}

    if ptm is not None:
        conf["ptm"] = float(_to_float_list(ptm)[0])
    if iptm is not None:
        conf["iptm"] = float(_to_float_list(iptm)[0])

    if pae is not None:
        arr = pae.detach().cpu().numpy() if hasattr(pae, "detach") else np.asarray(pae)
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 2:
            conf["pae"] = arr.tolist()

    return conf


def main():
    ap = argparse.ArgumentParser(description="Run one seeded ESMFold2 complex fold.")
    ap.add_argument("--input-json", required=True,
                    help="JSON with receptor/effector id+sequence.")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--repo", default="biohub/ESMFold2",
                    help="HF repo for the ESMFold2 head (matches esmfold2.def).")
    ap.add_argument("--num-loops", type=int, default=20)
    ap.add_argument("--num-sampling-steps", type=int, default=100)
    args = ap.parse_args()

    with open(args.input_json) as f:
        spec = json.load(f)
    rec, eff = spec["receptor"], spec["effector"]

    # Imports deferred so --help / unit tests don't need the GPU stack.
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder,
        ProteinInput,
        StructurePredictionInput,
    )
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

    model = ESMFold2Model.from_pretrained(args.repo).cuda().eval()
    spi = StructurePredictionInput(sequences=[
        ProteinInput(id=rec["id"], sequence=rec["sequence"]),
        ProteinInput(id=eff["id"], sequence=eff["sequence"]),
    ])
    result = ESMFold2InputBuilder().fold(
        model, spi,
        num_loops=args.num_loops,
        num_sampling_steps=args.num_sampling_steps,
        num_diffusion_samples=1,
        seed=args.seed,
    )

    os.makedirs(args.out_dir, exist_ok=True)
    cif_path = os.path.join(args.out_dir, "esmfold2_pred.cif")
    with open(cif_path, "w") as f:
        f.write(result.complex.to_mmcif())

    conf = build_confidences(
        result.plddt,
        ptm=getattr(result, "ptm", None),
        iptm=getattr(result, "iptm", getattr(result, "i_ptm", None)),
        pae=getattr(result, "pae", None),
    )
    with open(os.path.join(args.out_dir, "confidences.json"), "w") as f:
        json.dump(conf, f)

    if os.path.getsize(cif_path) < 500:
        print(f"ERROR: seed {args.seed} produced a suspiciously small mmCIF.",
              file=sys.stderr)
        sys.exit(1)

    import numpy as np
    print(f"seed {args.seed}: wrote {cif_path} "
          f"(plddt_mean={np.mean(conf['plddt']):.1f}, "
          f"ptm={conf.get('ptm')}, iptm={conf.get('iptm')})")


if __name__ == "__main__":
    main()
