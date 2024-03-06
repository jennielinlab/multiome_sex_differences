#!/usr/bin/env bash

workingdir="/projects/b1042/LinLab/realign-mismatch"

cat ${workingdir}/data/celltypes.txt | while read line
do

celltype=${line}

sbatch --job-name=DAR_${celltype} \
--error=${workingdir}/scripts/logs/${celltype}-DARs.err \
--output=${workingdir}/scripts/logs/${celltype}-DARs.out \
${workingdir}/scripts/runDARs.slurm ${celltype}

done
