#!/bin/bash
set -euo pipefail
input="data/XTal"
output_dir="data/PostProcessInputXTal"
smiles_dir="data/smiles"
lcse_dir="data/LCSEXTal"

python src/XTalinput.py --input_dir "$input" --output_dir "$output_dir" --smiles_dir "$smiles_dir" --LCSE_dir "$lcse_dir"

for folder in "$output_dir"/*;do
    pymol -cq -r src/cap_protein.py -- "$folder/protein.pdb" "$folder/protein_cap.pdb"
done
for folder in "$output_dir"/*;do
    python src/fixprotein.py "$folder/protein_cap.pdb" "$folder/protein_H_unordered.pdb"
done

for folder in "$output_dir"/*;do
    python src/reorder_pdb_residues.py "$folder/protein_H_unordered.pdb" "$folder/protein_H.pdb"
done

# Don't forget to run the verifyopenmm make it work