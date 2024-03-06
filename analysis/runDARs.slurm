#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=matthew.wrightpalmer@northwestern.edu
#SBATCH --mail-type=END,FAIL

module purge
module load R/4.1.1

celltype=${1}

workingdir="/projects/b1042/LinLab/realign-mismatch"

Rscript ${workingdir}/scripts/getDARs.R ${celltype}
