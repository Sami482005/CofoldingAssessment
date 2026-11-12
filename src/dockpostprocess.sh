#!/bin/bash
#!/bin/bash
#SBATCH --account=rrg-schapira
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=10:00:00
#SBATCH --job-name=AF3FF19
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --export=ALL
set -euo pipefail


# IF you face errors regarding to path or missing packages... just remove these 2 lines and uncomment the module load lines below.

export PYTHONNOUSERSITE=1
unset PYTHONPATH

module --force purge

if [ "${CONDA_DEFAULT_ENV:-}" != "dock-postprocess" ]; then
    echo "Expected the dock-postprocess Conda environment to be exported by sbatch." >&2
    echo "Submit this job with: sbatch --export=ALL src/dockpostprocess.sh" >&2
    exit 1
fi

if ! command -v dock-minimize >/dev/null 2>&1; then
    echo "dock-minimize is not available in the exported Conda environment." >&2
    echo "CONDA_PREFIX=${CONDA_PREFIX:-<unset>}" >&2
    exit 1
fi

source ~/miniconda3/etc/profile.d/conda.sh
conda activate dock-postprocess

path_to_work="/home/sami05/projects/rrg-schapira/sami05/cofoldingAssessment/data"

for input_dir in "$path_to_work/PostProcessInput"/*; do
    # Ensure it's a directory
    [ -d "$input_dir" ] || continue

    # Extract folder name (e.g., ICM_Object_10)
    obj_name=$(basename "$input_dir")

    # Define dedicated output directory for this object
    output_dir="$path_to_work/analysis/$obj_name"
    mkdir -p "$output_dir"

    echo "=== Processing $obj_name ==="

    # Run docking pipeline steps
    dock-init --input "$input_dir" --output "$output_dir" --receptor "$input_dir/protein_H.pdb"
    
    dock-minimize --results "$output_dir" --protein-forcefield ff19sb --solvent explicit --water-model opc --box-shape dodecahedron --platform CUDA
    dock-qc --results "$output_dir"
    dock-enrich --results "$output_dir"
    dock-landscape --results "$output_dir"
    dock-strain --results "$output_dir"
    dock-master --results "$output_dir"

    echo "Finished processing $obj_name"
    echo "----------------------------------------"
done