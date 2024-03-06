#!/usr/bin/env bash

workingdir="/projects/b1042/LinLab/NTS-bulk"

cat ${workingdir}/data/fastq_names.txt | while read line
do

sample=${line}

sbatch --job-name=trimAlign_${sample} \
--error=${workingdir}/scripts/logs/trim-align-${sample}.err \
--output=${workingdir}/scripts/logs/trim-align-${sample}.out \
${workingdir}/scripts/trim-and-align-single.slurm ${sample}

done
