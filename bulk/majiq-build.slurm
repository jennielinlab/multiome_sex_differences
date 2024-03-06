#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=matthew.wrightpalmer@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=majiq
#SBATCH --output=/projects/b1042/LinLab/NTS-bulk/scripts/logs/majiq-build-trimmed.out
#SBATCH --error=/projects/b1042/LinLab/NTS-bulk/scripts/logs/majiq-build-trimmed.err

cd /projects/b1042/LinLab/NTS-bulk

module purge
module load htslib
module load anaconda3
conda activate majiq-py3.8
source /home/mlw0278/condamajiq/bin/activate

export HTSLIB_LIBRARY_DIR=/software/supplemental/htslib-1.8/lib
export HTSLIB_INCLUDE_DIR=/software/supplemental/htslib-1.8/include

majiq build /home/mlw0278/reference/GencodeM23/gencode.vM23.primary_assembly.annotation.gff3 \
-c data/config-trimmed.txt -j 20 --min-intronic-cov 1 --simplify \
-o /projects/b1042/LinLab/NTS-bulk/output/majiq/trimmed/build/
  