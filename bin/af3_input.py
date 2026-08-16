#!/usr/bin/env python3
"""Generate AlphaFold 3 JSON input.

full   runs the data pipeline for both chains.
nomsa  uses --norun_data_pipeline with empty MSA fields.

Both modes are template-free, so no structural information from the reference
reaches the model. AF3's own PDB template search is disabled in the Nextflow
module by --max_template_date.

Usage:
    af3_input.py <receptor_seq> <effector_seq> <full|nomsa> [--output input.json]
"""

import argparse
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receptor_seq")
    parser.add_argument("effector_seq")
    parser.add_argument("mode", choices=["full", "nomsa"])
    parser.add_argument("--output", default="input.json")
    args = parser.parse_args()

    template_entries = []

    if args.mode == "full":
        # Full data pipeline (MSA + PDB template search) for BOTH chains.
        # No custom template injection: AF3 rejects a chain that has a custom
        # template alongside an auto-built MSA, so full mode is vanilla AF3.
        receptor_protein = {
            "id": "A",
            "sequence": args.receptor_seq,
        }
        effector_protein = {
            "id": "B",
            "sequence": args.effector_seq,
        }
        template_entries = []  # full mode never uses a custom template

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
