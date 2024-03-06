#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=36:00:00
#SBATCH --mail-user=matthew.wrightpalmer@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=fCounts
#SBATCH --output=/projects/b1042/LinLab/NTS-bulk/scripts/logs/featureCounts.out
#SBATCH --error=/projects/b1042/LinLab/NTS-bulk/scripts/logs/featureCounts.err

module purge
module load anaconda3
source activate featureCounts

workingdir="/projects/b1042/LinLab/NTS-bulk"

gtf="/home/mlw0278/reference/GencodeM23/gencode.vM23.primary_assembly.annotation.gtf"

allbams=$(ls -l ${workingdir}/output/aligned/trimmed/*.bam | awk '{print $9}' | tr "\n" " ")

outdir="${workingdir}/output/counts"
[[ -d ${outdir} ]] || mkdir ${outdir}

### Run featureCounts
### -p = paired end
### -T = 20 threads
### -s 2 = reverse stranded library
call="featureCounts -p -T 20 -s 2 -a ${gtf} \
-o ${outdir}/all-samples-counts-trimmed.txt \
${allbams}"

echo ${call}
eval ${call}
