#!/usr/bin/env python3
"""
extract_sequences.py
--------------------
Extract protein chain sequences from a PDB or mmCIF file.

Usage:
    python extract_sequences.py INPUT_PDB [--output sequences.json]

Output JSON:
    {
        "chains": [
            {"id": "A", "sequence": "MKFL...", "length": 245},
            {"id": "B", "sequence": "GTAL...", "length": 112}
        ],
        "chain_ids": ["A", "B"],
        "multimer_fasta_header": ">A_and_B",
        "multimer_fasta_sequence": "MKFL...:GTAL..."
    }

Also writes shell-sourceable variables to sequences.env:
    CHAIN_A_SEQ="MKFL..."
    CHAIN_B_SEQ="GTAL..."
    CHAIN_A_LEN=245
    CHAIN_B_LEN=112
    ALL_CHAIN_IDS="A B"
    NUM_CHAINS=2
"""

import argparse
import json
import os
import sys


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
    # Non-standard → X
    "MSE": "M", "HYP": "P", "TPO": "T", "SEP": "S", "PTR": "Y",
}


def extract_sequences_gemmi(path):
    """Extract sequences using gemmi (preferred — handles PDB and mmCIF)."""
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()

    chains = []
    seen = set()
    for model in st:
        for chain in model:
            if chain.name in seen:
                continue
            seq = []
            for res in chain.get_polymer():
                one = THREE_TO_ONE.get(res.name, "X")
                seq.append(one)
            if seq:
                seen.add(chain.name)
                chains.append({
                    "id": chain.name,
                    "sequence": "".join(seq),
                    "length": len(seq),
                })
        break  # first model only
    return chains


def extract_sequences_biopython(path):
    """Fallback: extract sequences using BioPython."""
    from Bio.PDB import PDBParser, MMCIFParser

    if path.endswith(".cif") or path.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("s", path)
    chains = []
    seen = set()
    for model in structure:
        for chain in model:
            if chain.id in seen:
                continue
            seq = []
            for res in chain:
                if res.id[0] != " ":
                    continue  # skip hetero
                one = THREE_TO_ONE.get(res.resname.strip(), "X")
                seq.append(one)
            if seq:
                seen.add(chain.id)
                chains.append({
                    "id": chain.id,
                    "sequence": "".join(seq),
                    "length": len(seq),
                })
        break
    return chains


def extract_sequences(path):
    """Try gemmi first, fall back to BioPython."""
    try:
        chains = extract_sequences_gemmi(path)
        if chains:
            return chains
    except ImportError:
        pass
    except Exception as e:
        print(f"WARNING: gemmi failed ({e}), trying BioPython...", file=sys.stderr)

    try:
        chains = extract_sequences_biopython(path)
        if chains:
            return chains
    except ImportError:
        pass
    except Exception as e:
        print(f"WARNING: BioPython failed ({e})", file=sys.stderr)

    print("ERROR: Neither gemmi nor BioPython could parse the file.", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Extract chain sequences from PDB/mmCIF")
    parser.add_argument("input", help="Input PDB or mmCIF file")
    parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
    parser.add_argument("--env", default=None, help="Output shell env file")
    parser.add_argument("--chains", nargs="+", default=None,
                        help="Chain IDs to use as CHAIN_A/CHAIN_B aliases "
                             "(e.g. --chains B C). Defaults to first two chains.")
    args = parser.parse_args()

    # If --chains not given on command line, fall back to BENCHMARK_CHAINS env var
    if not args.chains and os.environ.get("BENCHMARK_CHAINS"):
        args.chains = os.environ["BENCHMARK_CHAINS"].replace(",", " ").replace(":", " ").split()

    if not os.path.exists(args.input):
        print(f"ERROR: File not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    chains = extract_sequences(args.input)

    if not chains:
        print("ERROR: No protein chains found.", file=sys.stderr)
        sys.exit(1)

    result = {
        "chains": chains,
        "chain_ids": [c["id"] for c in chains],
        "num_chains": len(chains),
    }

    json_str = json.dumps(result, indent=2)

    if args.output:
        with open(args.output, "w") as f:
            f.write(json_str + "\n")
        print(f"Wrote {args.output} ({len(chains)} chains)", file=sys.stderr)
    else:
        print(json_str)

    # Write shell env file
    env_path = args.env or (args.output.replace(".json", ".env") if args.output else None)
    if env_path:
        with open(env_path, "w") as f:
            f.write("# Auto-generated by extract_sequences.py\n")
            f.write(f'NUM_CHAINS={len(chains)}\n')
            f.write(f'ALL_CHAIN_IDS="{" ".join(c["id"] for c in chains)}"\n')
            for c in chains:
                f.write(f'CHAIN_{c["id"]}_SEQ="{c["sequence"]}"\n')
                f.write(f'CHAIN_{c["id"]}_LEN={c["length"]}\n')
            # Convenience aliases: use --chains if specified, else first two
            alias_chains = []
            if args.chains:
                # Accept both space-separated ("B" "C") and comma-separated ("B,C")
                requested = []
                for item in args.chains:
                    requested.extend(item.split(","))
                requested = [c.strip() for c in requested if c.strip()]
                chain_map = {c["id"]: c for c in chains}
                for cid in requested[:2]:
                    if cid in chain_map:
                        alias_chains.append(chain_map[cid])
                    else:
                        print(f"WARNING: --chains {cid} not found in PDB (available: {[c['id'] for c in chains]})", file=sys.stderr)
            if not alias_chains and len(chains) >= 2:
                alias_chains = chains[:2]
            if len(alias_chains) >= 2:
                f.write(f'\n# Convenience aliases (chains {alias_chains[0]["id"]} and {alias_chains[1]["id"]})\n')
                f.write(f'CHAIN_A_SEQ="{alias_chains[0]["sequence"]}"\n')
                f.write(f'CHAIN_B_SEQ="{alias_chains[1]["sequence"]}"\n')
                f.write(f'CHAIN_A_LEN={alias_chains[0]["length"]}\n')
                f.write(f'CHAIN_B_LEN={alias_chains[1]["length"]}\n')
        print(f"Wrote {env_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
