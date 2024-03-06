#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=48:00:00
#SBATCH --mail-user=matthew.wrightpalmer@northwestern.edu
#SBATCH --mail-type=END,FAIL

module purge
module load anaconda3
source activate cutadapt

module load STAR
module load samtools

workingdir="/projects/b1042/LinLab/NTS-bulk"
trimdir="${workingdir}/output/trimmed"

indexDir="/projects/b1158/star-indexes/GencodeM23"

aligndir="${workingdir}/output/aligned/trimmed"
[[ -d ${aligndir} ]] || mkdir ${aligndir}

sample=${1}

fqdir="${workingdir}/data"

call="cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
-A ACTGTCTCTTATACACATCTGACGCTGCCGACGA \
-j 20 --trim-n -m 25 -q 10,10 --nextseq-trim=10 \
-o ${trimdir}/${sample}_R1_001.fastq.gz \
-p ${trimdir}/${sample}_R2_001.fastq.gz \
${fqdir}/${sample}_R1_001.fastq.gz \
${fqdir}/${sample}_R2_001.fastq.gz"

echo ${call}
eval ${call} 2>&1 | tee ${trimdir}/logs/${sample}.log

### Run STAR
### Two pass mode: --twopassMode Basic \
### Flags recommended for cufflinks/stringtie and used in Gitler:
# --outFilterIntronMotifs RemoveNoncanonical \
# --alignEndsType EndToEnd
### ENCODE options from STAR manual:
# --outFilterType BySJout
# --outFilterMultimapNmax 20
# --alignSJoverhangMin 8
# --alignSJDBoverhangMin 1
# --outFilterMismatchNmax 999
# --outFilterMismatchNoverReadLmax 0.04
# --alignIntronMin 20
# --alignIntronMax 1000000
# --alignMatesGapMax 1000000

call="STAR --runThreadN 20 \
--genomeDir ${indexDir} \
--readFilesIn ${trimdir}/${sample}_R1_001.fastq.gz ${trimdir}/${sample}_R2_001.fastq.gz \
--readFilesCommand zcat \
--twopassMode Basic \
--outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--alignEndsType EndToEnd \
--quantMode GeneCounts \
--outFileNamePrefix ${aligndir}/${sample}. \
--outSAMattributes NH HI AS nM NM"

echo ${call}
eval ${call}

#rm ${tmpdir}/*.fq.gz
#rm ${tmptrimdir}/*.fq.gz


samtools index -@ 16 ${aligndir}/${sample}.Aligned.sortedByCoord.out.bam
