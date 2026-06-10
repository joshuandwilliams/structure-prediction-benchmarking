#!/usr/bin/env python3
"""
extract_sequences_from_fasta.py
-------------------------------
Alternative input adapter for the structure-prediction benchmark pipeline.

Reads a FASTA file containing exactly two protein entries (receptor + effector)
and emits a sequences.json file in the SAME shape that extract_sequences.py
produces from a reference PDB, so all downstream Nextflow processes remain
unchanged.

Usage:
    python extract_sequences_from_fasta.py INPUT_FASTA \
        --output sequences.json \
        --env sequences.env \
        --chains A B

The two --chains values become the chain IDs assigned to the receptor and
effector respectively (defaults: A B).  They are purely synthetic labels
used by Boltz/Chai/etc. as chain identifiers — there is no reference PDB
in FASTA mode, so no real PDB chain needs to match.

Output JSON (matches extract_sequences.py):
    {
        "chains": [
            {"id": "A", "sequence": "MKFL...", "length": 245},
            {"id": "B", "sequence": "GTAL...", "length": 112}
        ],
        "chain_ids": ["A", "B"],
        "num_chains": 2
    }

Also writes shell-sourceable variables to sequences.env, with the same
CHAIN_*_SEQ / CHAIN_A_SEQ alias scheme as extract_sequences.py.
"""

import argparse
import json
import os
import re
import sys

# Standard 20 + Sec/Pyl + ambiguity codes accepted in protein FASTAs.
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYUOBZXJ*-")


def parse_fasta(path):
    """
    Minimal FASTA parser.  Returns a list of (header, sequence) tuples in
    file order.  Sequence whitespace is stripped; case is upper-cased;
    trailing '*' stop codons are removed.  Raises ValueError on malformed
    input.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"FASTA file not found: {path}")

    entries = []
    current_header = None
    current_seq_chunks = []

    with open(path) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    entries.append((current_header, "".join(current_seq_chunks)))
                current_header = line[1:].strip() or f"entry_{len(entries) + 1}"
                current_seq_chunks = []
            else:
                if current_header is None:
                    raise ValueError(
                        f"FASTA parse error at line {lineno}: sequence data "
                        f"appears before any '>' header."
                    )
                # Strip any internal whitespace, upper-case.
                cleaned = re.sub(r"\s+", "", line).upper()
                current_seq_chunks.append(cleaned)

    if current_header is not None:
        entries.append((current_header, "".join(current_seq_chunks)))

    # Drop trailing stop codons that some tools embed.
    cleaned_entries = []
    for header, seq in entries:
        seq = seq.rstrip("*")
        cleaned_entries.append((header, seq))

    return cleaned_entries


def validate_sequence(header, seq):
    """Sanity-check that a sequence is non-empty and contains only AA chars."""
    if not seq:
        raise ValueError(f"Entry '{header}' has empty sequence.")
    bad = set(seq) - VALID_AA
    if bad:
        raise ValueError(
            f"Entry '{header}' contains non-amino-acid characters: "
            f"{sorted(bad)}.  Are you sure this is a protein FASTA?"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Extract two protein chain sequences from a FASTA file "
                    "(receptor + effector) for the benchmark pipeline."
    )
    parser.add_argument("input", help="Input FASTA file with exactly two entries")
    parser.add_argument("--output", default=None,
                        help="Output JSON file (default: stdout)")
    parser.add_argument("--env", default=None, help="Output shell env file")
    parser.add_argument(
        "--chains", nargs="+", default=None,
        help="Chain IDs to assign to the two FASTA entries, in order "
             "(receptor effector).  Defaults to 'A B'.  These are synthetic "
             "labels used by the predictors as chain identifiers.",
    )
    args = parser.parse_args()

    # Same env-var fallback convention as extract_sequences.py.
    if not args.chains and os.environ.get("BENCHMARK_CHAINS"):
        args.chains = (
            os.environ["BENCHMARK_CHAINS"]
            .replace(",", " ").replace(":", " ").split()
        )
    if not args.chains:
        args.chains = ["A", "B"]

    # Normalise --chains: accept both "A B" and "A,B" forms.
    requested_ids = []
    for item in args.chains:
        requested_ids.extend(item.split(","))
    requested_ids = [c.strip() for c in requested_ids if c.strip()]
    if len(requested_ids) < 2:
        print(
            f"ERROR: Need exactly two chain IDs for --chains, got "
            f"{requested_ids}.",
            file=sys.stderr,
        )
        sys.exit(1)
    receptor_id, effector_id = requested_ids[0], requested_ids[1]
    if receptor_id == effector_id:
        print(
            f"ERROR: Receptor and effector chain IDs must differ "
            f"(both were '{receptor_id}').",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        entries = parse_fasta(args.input)
    except (FileNotFoundError, ValueError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if len(entries) != 2:
        print(
            f"ERROR: FASTA must contain exactly 2 entries (receptor + "
            f"effector), found {len(entries)} in {args.input}.",
            file=sys.stderr,
        )
        for i, (h, s) in enumerate(entries, 1):
            print(f"  entry {i}: >{h}  ({len(s)} aa)", file=sys.stderr)
        sys.exit(1)

    (rec_header, rec_seq), (eff_header, eff_seq) = entries
    try:
        validate_sequence(rec_header, rec_seq)
        validate_sequence(eff_header, eff_seq)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    chains = [
        {"id": receptor_id, "sequence": rec_seq, "length": len(rec_seq),
         "fasta_header": rec_header},
        {"id": effector_id, "sequence": eff_seq, "length": len(eff_seq),
         "fasta_header": eff_header},
    ]

    result = {
        "chains": chains,
        "chain_ids": [c["id"] for c in chains],
        "num_chains": len(chains),
        "source": "fasta",
        "source_file": os.path.basename(args.input),
    }

    json_str = json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w") as f:
            f.write(json_str + "\n")
        print(
            f"Wrote {args.output} (2 chains from FASTA: "
            f"{receptor_id}={len(rec_seq)}aa, {effector_id}={len(eff_seq)}aa)",
            file=sys.stderr,
        )
    else:
        print(json_str)

    # Shell env file — same convention as extract_sequences.py so anything
    # downstream that sources sequences.env still works.
    env_path = args.env or (
        args.output.replace(".json", ".env") if args.output else None
    )
    if env_path:
        with open(env_path, "w") as f:
            f.write("# Auto-generated by extract_sequences_from_fasta.py\n")
            f.write(f"NUM_CHAINS={len(chains)}\n")
            f.write(f'ALL_CHAIN_IDS="{" ".join(c["id"] for c in chains)}"\n')
            for c in chains:
                f.write(f'CHAIN_{c["id"]}_SEQ="{c["sequence"]}"\n')
                f.write(f'CHAIN_{c["id"]}_LEN={c["length"]}\n')
            f.write(
                f"\n# Convenience aliases (chains {receptor_id} and "
                f"{effector_id})\n"
            )
            f.write(f'CHAIN_A_SEQ="{rec_seq}"\n')
            f.write(f'CHAIN_B_SEQ="{eff_seq}"\n')
            f.write(f"CHAIN_A_LEN={len(rec_seq)}\n")
            f.write(f"CHAIN_B_LEN={len(eff_seq)}\n")
        print(f"Wrote {env_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
