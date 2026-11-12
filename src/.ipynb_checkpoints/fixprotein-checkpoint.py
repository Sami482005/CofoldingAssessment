import argparse
from collections import OrderedDict
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from pathlib import Path

ATOM_RENAMES = {
    "ACE": {
        "1HH3": "H1",
        "2HH3": "H2",
        "3HH3": "H3",
        "HH31": "H1",
        "HH32": "H2",
        "HH33": "H3",
    },
    "NME": {
        "1HH3": "H1",
        "2HH3": "H2",
        "3HH3": "H3",
        "HH31": "H1",
        "HH32": "H2",
        "HH33": "H3",
        "CH3": "C",
    },
}

TERMINAL_RESIDUES = {"NME", "UNL"}


def rename_atom(line):
    """Rename terminal atoms while preserving PDB fixed-column formatting."""
    residue_name = line[17:21].strip()
    atom_name = line[12:16].strip()

    replacement = ATOM_RENAMES.get(residue_name, {}).get(atom_name)
    if replacement is None:
        return line

    # PDB atom names occupy columns 13-16 and must remain four characters wide.
    replacement = replacement.ljust(4)
    return line[:12] + replacement + line[16:]


def process_pdb(input_path, output_path):
    """
    Rename ACE/NME atoms and insert TER records after terminal residues.

    Only ATOM and HETATM records are processed. Other records are preserved.
    """
    residues = OrderedDict()
    residue_order = []

    with open(input_path, "r") as input_file:
        for line in input_file:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            residue_name = line[17:21].strip()
            chain_id = line[21]
            residue_number = line[22:26]
            insertion_code = line[26]

            residue_key = (
                chain_id,
                residue_number,
                insertion_code,
            )

            if residue_key not in residues:
                residues[residue_key] = {
                    "residue_name": residue_name,
                    "lines": [],
                }
                residue_order.append(residue_key)

            residues[residue_key]["lines"].append(rename_atom(line))

    output_lines = []

    for residue_key in residue_order:
        residue = residues[residue_key]

        output_lines.extend(residue["lines"])

        if residue["residue_name"] in TERMINAL_RESIDUES:
            output_lines.append("TER\n")

    with open(output_path, "w") as output_file:
        output_file.writelines(output_lines)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Rename ACE/NME terminal atoms and add TER records after "
            "NME or UNL residues."
        )
    )
    parser.add_argument("input_pdb", help="Input PDB file")
    parser.add_argument("output_pdb", help="Output PDB file")

    args = parser.parse_args()
    path = prep_protein(args.input_pdb)

    process_pdb(
        input_path=path,
        output_path=args.output_pdb,
    )

def prep_protein(protein_pdb_path):
    protein_pdb_path = Path(protein_pdb_path)

    fixer = PDBFixer(filename=str(protein_pdb_path))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.addMissingHydrogens(7.0)

    output_protein_pdb_path = protein_pdb_path.parent / "protein_H.pdb"

    with open(output_protein_pdb_path, "w") as f:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True,
        )

    print(f"Protein written: {output_protein_pdb_path}")
    return output_protein_pdb_path

if __name__ == "__main__":
    main()