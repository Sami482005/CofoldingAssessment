#!/usr/bin/env python3
"""
reorder_pdb_residues.py
=======================

Reorder the residues of a PDB file so they appear in proper sequence order.

Typical use case: capping / protonation tools (e.g. adding ACE/NME caps) often
append or prepend the new cap residues at the top or bottom of the file, so the
residue order no longer matches the numbering, e.g.::

    ACE 950, NME 1112, ACE 1114, SER 951, GLY 952, ...   (wrong)
    ACE 950, SER 951, ..., SER 1111, NME 1112, TER, ACE 1114, MET 1115, ...

What the script does
--------------------
1. Groups ATOM/HETATM records into residues keyed by
   (chain ID, residue number, insertion code).
2. Keeps chains in the order they first appear in the file.
3. Sorts residues within each chain by (residue number, insertion code).
   Atom order *inside* each residue is preserved exactly.
4. Renumbers atom serial numbers from 1 and remaps any CONECT records so
   bonds still point at the correct atoms.
5. Writes a TER record after each chain and, optionally (``--ter-at-gaps``),
   wherever the residue numbering jumps (a chain break, e.g. between the
   NME cap at 1112 and the ACE cap at 1114). This stops downstream tools
   from bonding across the gap.

Header records (CRYST1, REMARK, ...) that appear before the first atom are
kept at the top; old TER/END/MODEL records are dropped and regenerated.
Multi-model files (NMR ensembles) are not supported.

Usage
-----
    python reorder_pdb_residues.py input.pdb output.pdb [--ter-at-gaps]
"""

import argparse
import sys
from collections import OrderedDict

ATOM_RECORDS = ("ATOM  ", "HETATM")
DROPPED_RECORDS = ("TER", "END", "MODEL", "ENDMDL", "CONECT", "MASTER")


def residue_key(line):
    """Return (chain, resseq, icode) for an ATOM/HETATM line."""
    chain = line[21]
    resseq = int(line[22:26])
    icode = line[26]
    return chain, resseq, icode


def parse_pdb(path):
    """
    Read a PDB file.

    Returns
    -------
    header : list[str]   non-atom records that precede the first atom
    chains : OrderedDict[chain -> OrderedDict[(resseq, icode) -> list[str]]]
    conect : list[str]   original CONECT lines
    """
    header, conect = [], []
    chains = OrderedDict()
    seen_atom = False

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            record = line[:6]

            if record in ATOM_RECORDS:
                seen_atom = True
                chain, resseq, icode = residue_key(line)
                chains.setdefault(chain, OrderedDict()) \
                      .setdefault((resseq, icode), []).append(line)
            elif record.startswith("CONECT"):
                conect.append(line)
            elif record.strip().startswith(DROPPED_RECORDS):
                if record.startswith("ENDMDL"):
                    sys.exit("Error: multi-model PDB files are not supported.")
                continue
            elif not seen_atom:
                header.append(line)
            # Any other record appearing between atoms is dropped.

    if not chains:
        sys.exit(f"Error: no ATOM/HETATM records found in {path}")
    return header, chains, conect


def ter_line(serial, last_atom_line):
    """Build a TER record that follows the given atom line."""
    return f"TER   {serial:>5d}      {last_atom_line[17:27]}"


def write_sorted(header, chains, conect, out_path, ter_at_gaps):
    """Write residues in sorted order with renumbered serials and TER records."""
    out = list(header)
    serial_map = {}   # old atom serial -> new atom serial
    serial = 0

    for chain_residues in chains.values():
        sorted_keys = sorted(chain_residues)  # (resseq, icode)
        prev_resseq = None
        last_atom = None

        for resseq, icode in sorted_keys:
            # Optional TER at a numbering gap (chain break).
            if ter_at_gaps and prev_resseq is not None and resseq > prev_resseq + 1:
                serial += 1
                out.append(ter_line(serial, last_atom))

            for atom in chain_residues[(resseq, icode)]:
                serial += 1
                old = atom[6:11].strip()
                if old:
                    serial_map[old] = serial
                atom = f"{atom[:6]}{serial:>5d}{atom[11:]}"
                out.append(atom)
                last_atom = atom
            prev_resseq = resseq

        # TER at the end of every chain.
        serial += 1
        out.append(ter_line(serial, last_atom))

    # Remap CONECT records (columns: 7-11, 12-16, 17-21, ...).
    for line in conect:
        fields = [line[i:i + 5].strip() for i in range(6, len(line), 5)]
        fields = [f for f in fields if f]
        try:
            new = [serial_map[f] for f in fields]
        except KeyError:
            continue  # bond to an atom that no longer exists; skip it
        out.append("CONECT" + "".join(f"{n:>5d}" for n in new))

    out.append("END")
    with open(out_path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Reorder PDB residues by chain and residue number.")
    parser.add_argument("input", help="input PDB file")
    parser.add_argument("output", help="output PDB file")
    parser.add_argument("--ter-at-gaps", action="store_true",
                        help="insert a TER record wherever residue numbering "
                             "jumps (chain break)")
    args = parser.parse_args()

    header, chains, conect = parse_pdb(args.input)
    write_sorted(header, chains, conect, args.output, args.ter_at_gaps)
    n_res = sum(len(r) for r in chains.values())
    print(f"Wrote {args.output}: {len(chains)} chain(s), {n_res} residues.")


if __name__ == "__main__":
    main()
