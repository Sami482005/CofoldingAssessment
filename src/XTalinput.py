import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import os



def write_protein_only_pdb(input_pdb: Path, output_pdb: Path) -> None:
    """Write only ATOM records from the input PDB file."""
    with input_pdb.open("r") as input_file, output_pdb.open("w") as output_file:
        for line in input_file:
            if line.startswith("ATOM"):
                output_file.write(line)

        output_file.write("END\n")


def write_ligand_only_pdb(input_pdb: Path, output_pdb: Path) -> None:
    """Write HETATM records from the input PDB file."""
    with input_pdb.open("r") as input_file, output_pdb.open("w") as output_file:
        for line in input_file:
            if line.startswith("HETATM"):
                output_file.write(line)

        output_file.write("END\n")


def write_outputs(
    input_dir: str,
    output_dir: str,
    smiles_dir: str,
    LCSE_dir: str
) -> None:
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    smiles_path = Path(smiles_dir)
    lcse_path = Path(LCSE_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    lcse_path.mkdir(parents=True, exist_ok=True)

    for pdb_file in sorted(input_path.glob("*.pdb")):
        molecule_name = pdb_file.stem
        molecule_output_dir = output_path / molecule_name
        molecule_output_dir.mkdir(parents=True, exist_ok=True)
        write_protein_only_pdb(pdb_file, molecule_output_dir / "protein.pdb")
        write_ligand_only_pdb(pdb_file, molecule_output_dir / "ligand.pdb")
        prep_ligand(molecule_output_dir, smiles_path, lcse_path, molecule_name)


def assign_bond_orders(template_smiles: str, ligand: Chem.Mol, molecule_name: str) -> Chem.Mol:
    """Transfer bond orders and formal charges from a SMILES template onto a PDB ligand.

    First tries a match that includes explicit hydrogens, which keeps the
    crystal hydrogen coordinates. If the two protonation states disagree
    (e.g. the PDB protonates an atom that the SMILES leaves charged), that
    match fails. In that case the match is retried on heavy atoms only, and
    the hydrogens are later rebuilt from the SMILES protonation state.
    """
    template = Chem.MolFromSmiles(template_smiles)
    if template is None:
        raise ValueError(f"Could not parse SMILES for {molecule_name}: {template_smiles}")

    # Attempt 1: match with hydrogens on both sides.
    try:
        return AllChem.AssignBondOrdersFromTemplate(Chem.AddHs(template), ligand)
    except ValueError:
        pass

    # Attempt 2: heavy atoms only, so differing protonation no longer matters.
    print(
        f"Hydrogens of {molecule_name} differ between PDB and SMILES; "
        "matching on heavy atoms and rebuilding hydrogens from the SMILES."
    )
    try:
        return AllChem.AssignBondOrdersFromTemplate(
            Chem.RemoveHs(template),
            Chem.RemoveHs(ligand, sanitize=False),
        )
    except ValueError as error:
        raise ValueError(
            f"Failed to assign bond orders for {molecule_name}: {error}"
        ) from error


def prep_ligand(molecule_output_dir: Path, smiles_path: Path, lcse_path: Path, molecule_name: str) -> None:
    ligand = Chem.MolFromPDBFile(str(molecule_output_dir / "ligand.pdb"), removeHs=False)
    if ligand is None:
        print(f"Failed to read ligand from {molecule_output_dir.name}")
        return

    smiles_file = smiles_path / f"{molecule_name}.txt"
    if not smiles_file.exists():
        raise FileNotFoundError(
            f"SMILES file not found: {smiles_file}"
        )

    smiles = smiles_file.read_text().strip()
    mol = assign_bond_orders(smiles, ligand, molecule_name)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)  # Add hydrogens to the molecule
    sdf_path = molecule_output_dir / "ligand.sdf"
    lcse_sdf_path = lcse_path / f"{molecule_name}.sdf"

    with Chem.SDWriter(str(sdf_path)) as writer:
        writer.write(mol)

    with Chem.SDWriter(str(lcse_sdf_path)) as writer:
        writer.write(mol)
    
    os.remove(molecule_output_dir / "ligand.pdb")  # Remove the temporary ligand PDB file
    print(f"Ligand written: {sdf_path}")
        

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare input files for Dock-PostProcess"
    )

    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing input PDB files",
    )

    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output folders will be written",
    )

    parser.add_argument(
        "--smiles_dir",
        required=True,
        help="Directory containing SMILES files for each molecule",
    )

    parser.add_argument(
        "--LCSE_dir",
        required=True,
        help="Directory where LCSE output files will be written",
    )

    args = parser.parse_args()

    write_outputs(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        smiles_dir=args.smiles_dir,
        LCSE_dir=args.LCSE_dir
    )

if __name__ == "__main__":
    main()

# python src/XTalinput.py --input_dir data/XTal --output_dir data/PostProcessInputXTal --smiles_dir data/smiles --LCSE_dir data/LCSEXTal
