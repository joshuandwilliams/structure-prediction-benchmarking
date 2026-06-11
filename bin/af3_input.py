#!/usr/bin/env python3
"""
Generate AlphaFold 3 JSON input.

Two modes:
  full   — let AF3 run its full data pipeline (MSA + PDB template search).
            An effector structural template is injected into the JSON for the
            effector chain (chain B), which causes AF3 to use it instead of
            searching the PDB for effector templates.  The receptor (chain A)
            still gets the full MSA + template treatment.
  nomsa  — --norun_data_pipeline mode: MSA fields are empty and no PDB
            template search is run.  An effector template is still injected
            if provided.

Usage:
    af3_input.py <receptor_seq> <effector_seq> <mode> \\
        [--template <effector_template.cif>] [--output input.json]
"""

import argparse
import json
import os
import sys


def load_template(mmcif_path, n_effector_residues):
    """Return a list with one AF3 template entry, or [] if path is empty/missing."""
    if not mmcif_path:
        return []
    if not os.path.exists(mmcif_path) or os.path.getsize(mmcif_path) == 0:
        return []
    mmcif_str = open(mmcif_path).read()
    # Empty queryIndices/templateIndices: AF3 performs its own sequence
    # alignment between the query and template, which is more robust than
    # an explicit identity mapping when residues are missing from the crystal.
    return [
        {
            "mmcif": mmcif_str,
            "queryIndices": [],
            "templateIndices": [],
        }
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receptor_seq")
    parser.add_argument("effector_seq")
    parser.add_argument("mode", choices=["full", "nomsa"])
    parser.add_argument("--template", default=None,
                        help="Path to effector mmCIF template (empty = no template)")
    parser.add_argument("--output", default="input.json")
    args = parser.parse_args()

    n_eff = len(args.effector_seq)
    template_entries = load_template(args.template, n_eff)

    if args.mode == "full":
        # Full data pipeline: MSA + PDB template search for receptor.
        # Effector gets the injected template (overrides PDB search for that chain).
        receptor_protein = {
            "id": "A",
            "sequence": args.receptor_seq,
        }
        effector_protein = {
            "id": "B",
            "sequence": args.effector_seq,
        }
        if template_entries:
            effector_protein["templates"] = template_entries

    else:  # nomsa
        # No data pipeline: empty MSA, no PDB template search.
        # Effector gets the injected template if provided.
        receptor_protein = {
            "id": "A",
            "sequence": args.receptor_seq,
            "unpairedMsa": "",
            "pairedMsa": "",
            "templates": [],
        }
        effector_protein = {
            "id": "B",
            "sequence": args.effector_seq,
            "unpairedMsa": "",
            "pairedMsa": "",
            "templates": template_entries,
        }

    data = {
        "name": "benchmark_test",
        "modelSeeds": [42, 123, 456, 789, 1024],
        "dialect": "alphafold3",
        "version": 1,
        "sequences": [
            {"protein": receptor_protein},
            {"protein": effector_protein},
        ],
    }

    with open(args.output, "w") as f:
        json.dump(data, f, indent=2)

    mode_note = f"{args.mode} mode, {'with' if template_entries else 'no'} effector template"
    print(f"Written AF3 JSON ({mode_note}): {args.output}")


if __name__ == "__main__":
    main()
