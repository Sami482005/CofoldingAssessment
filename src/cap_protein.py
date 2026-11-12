"""Cap protein termini using the same logic as prepare_protein.py.

This script is intended to be run with the PyMOL Python interpreter, for
example:

    pymol -cq cap_protein.py -- input.pdb capped.pdb

It can also be imported from a PyMOL Python session:

    from cap_protein import cap_protein
    cap_protein("input.pdb", "capped.pdb")
"""

from collections import defaultdict
import argparse
import logging
from pathlib import Path

import numpy as np
from pymol import cmd, editor


LOGGER = logging.getLogger(__name__)


def _chain_selection(chain):
    """Return the chain clause used by the original script."""
    return f" and chain {chain}" if chain else ""


def _find_segments(residue_numbers):
    """Find segment termini using the original residue-gap rule."""
    list_base = np.sort(list(set(residue_numbers)))
    list1 = np.concatenate([[0], list_base])
    list2 = np.concatenate([list_base, [1000000]])
    idxs = np.where(list2 > list1 + 1)

    n_ter = [list2[i] for i in idxs[0] if list2[i] < 1e5]
    c_ter = [list1[i] for i in idxs[0] if list1[i] > 0]

    if len(n_ter) != len(c_ter):
        if len(n_ter) < len(c_ter):
            n_ter = [list_base[0]] + n_ter
        else:
            c_ter = c_ter + [list_base[-1]]

    return n_ter, c_ter


def cap_protein(input_pdb, output_pdb, object_name="prot"):
    """Attach ACE and NME caps to all protein segments in a PDB file.

    Protein segments are split exactly as in the original script: a gap of
    more than six between sorted residue numbers starts a new segment.
    Existing NME handling is also preserved.
    """
    input_path = Path(input_pdb)
    output_path = Path(output_pdb)
    if not input_path.is_file():
        raise FileNotFoundError(f"Input PDB file does not exist: {input_path}")

    cmd.set("retain_order", 1)
    cmd.load(str(input_path), object=object_name)
    cmd.remove(f"{object_name} and hydro")

    chain_residues = defaultdict(list)
    stored = []
    cmd.iterate(
        f"{object_name} and polymer.protein",
        "stored.append((int(resi), chain))",
        space={"stored": stored},
    )
    if not stored:
        raise ValueError("The input PDB contains no polymer protein atoms.")

    for residue_number, chain in stored:
        chain_residues[chain].append(residue_number)

    for chain, residue_numbers in chain_residues.items():
        chain_query = _chain_selection(chain)
        n_ter, c_ter = _find_segments(residue_numbers)

        for residue_number in n_ter:
            selection = (
                f"{object_name} and resi {residue_number} and name N"
                f"{chain_query}"
            )
            LOGGER.info("Adding ACE at %s", selection)
            editor.attach_amino_acid(selection, "ace")

        for residue_number in c_ter:
            if cmd.select(
                "c_ter_check",
                f"{object_name} and resi {residue_number + 1} "
                f"and resn NME{chain_query}",
            ):
                LOGGER.info("Skipping already capped C-terminus at residue %s",
                            residue_number)
                continue
            if cmd.select(
                "c_ter_check",
                f"{object_name} and resi {residue_number} and name OXT"
                f"{chain_query}",
            ):
                LOGGER.info("Skipping real C-terminus at residue %s",
                            residue_number)
                continue
            if cmd.select(
                "c_ter_check",
                f"{object_name} and resi {residue_number} and name HXT"
                f"{chain_query}",
            ):
                LOGGER.info("Skipping real H-terminus at residue %s",
                            residue_number)
                continue

            selection = (
                f"{object_name} and resi {residue_number} and name C"
                f"{chain_query}"
            )
            LOGGER.info("Adding NME at %s", selection)
            editor.attach_amino_acid(selection, "nme")

    cmd.save(str(output_path), object_name)
    cmd.delete(object_name)


def main():
    parser = argparse.ArgumentParser(
        description="Attach ACE/NME caps to protein termini in a PDB file."
    )
    parser.add_argument("input_pdb", help="Input PDB containing protein atoms")
    parser.add_argument("output_pdb", help="Output PDB to write")
    parser.add_argument(
        "--object-name",
        default="prot",
        help="Temporary PyMOL object name (default: prot)",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    cap_protein(args.input_pdb, args.output_pdb, args.object_name)


if __name__ == "__main__":
    main()
