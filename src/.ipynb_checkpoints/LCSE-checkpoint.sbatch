#!/bin/bash
#!/bin/bash
#SBATCH --account=rrg-schapira
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --job-name=LCSECPCM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x.out
#SBATCH --export=ALL
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate lcse

# python lcse/LCSE.py --ligand-dir data/LCSEinput --output-dir data/LCSEOutput

python lcse/LCSE.py --ligand-dir data/LCSEXTal --output-dir data/LCSEOutput
