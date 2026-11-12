#!/bin/bash
#SBATCH --account=rrg-schapira
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --job-name=LCSECPCMXtalNew
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x_%j_Nov102026.out
#SBATCH --export=ALL
set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate PCCL

# python lcse/LCSE.py --ligand-dir data/LCSEinput --output-dir data/LCSEOutput

python /home/sami05/projects/rrg-schapira/sami05/lcse/LCSE.py --ligand-dir data/LCSEXTal --output-dir data/LCSEOutput
