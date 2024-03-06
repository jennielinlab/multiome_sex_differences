#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --mail-user=matthew.wrightpalmer@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=majiq
#SBATCH --output=/projects/b1042/LinLab/NTS-bulk/scripts/logs/majiq-deltapsi-trimmed.out
#SBATCH --error=/projects/b1042/LinLab/NTS-bulk/scripts/logs/majiq-deltapsi-trimmed.err


cd /projects/b1042/LinLab/NTS-bulk/output/majiq/trimmed/build

module purge
module load htslib
module load anaconda3
conda activate majiq-py3.8
source /home/mlw0278/condamajiq/bin/activate

export HTSLIB_LIBRARY_DIR=/software/supplemental/htslib-1.8/lib
export HTSLIB_INCLUDE_DIR=/software/supplemental/htslib-1.8/include

outputDir="/projects/b1042/LinLab/NTS-bulk/output/majiq/trimmed/deltapsi"

### Control Young, MvF
majiq deltapsi \
-grp1 Y-F-C-1_S1.Aligned.sortedByCoord.out.majiq \
Y-F-C-5_S2.Aligned.sortedByCoord.out.majiq \
Y-F-C--8_S3.Aligned.sortedByCoord.out.majiq \
Y-F-C-9_S4.Aligned.sortedByCoord.out.majiq \
-grp2 Y-M-C-25_S17.Aligned.sortedByCoord.out.majiq \
Y-M-C-27_S18.Aligned.sortedByCoord.out.majiq \
Y-M-C-31_S19.Aligned.sortedByCoord.out.majiq \
Y-M-C-33_S20.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FYC MYC

### NTS Young, MvF 
majiq deltapsi \
-grp1 Y-F-NTS-2_S5.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-4_S6.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-6_S7.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-7_S8.Aligned.sortedByCoord.out.majiq \
-grp2 Y-M-NTS-26_S21.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-28_S22.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-30_S23.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-34_S24.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FYN MYN

### Control Aged, MvF
majiq deltapsi \
-grp1 OLD-F-C-13_S9.Aligned.sortedByCoord.out.majiq \
OLD-F-C-19_S10.Aligned.sortedByCoord.out.majiq \
OLD-F-C-21_S11.Aligned.sortedByCoord.out.majiq \
OLD-F-C-23_S12.Aligned.sortedByCoord.out.majiq \
-grp2 OLD-M-C-37_S25.Aligned.sortedByCoord.out.majiq \
OLD-M-C-41_S26.Aligned.sortedByCoord.out.majiq \
OLD-M-C-43_S27.Aligned.sortedByCoord.out.majiq \
OLD-M-C-45_S28.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FAC MAC

### NTS Aged, MvF
majiq deltapsi \
-grp1 OLD-F-NTS-14_S13.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-18_S14.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-20_S15.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-24_S16.Aligned.sortedByCoord.out.majiq \
-grp2 Aged-M-NTS-38_S1.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-42_S2.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-44_S3.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-46_S4.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FAN MAN

### Female Control, young vs aged
majiq deltapsi \
-grp1 Y-F-C-1_S1.Aligned.sortedByCoord.out.majiq \
Y-F-C-5_S2.Aligned.sortedByCoord.out.majiq \
Y-F-C--8_S3.Aligned.sortedByCoord.out.majiq \
Y-F-C-9_S4.Aligned.sortedByCoord.out.majiq \
-grp2 OLD-F-C-13_S9.Aligned.sortedByCoord.out.majiq \
OLD-F-C-19_S10.Aligned.sortedByCoord.out.majiq \
OLD-F-C-21_S11.Aligned.sortedByCoord.out.majiq \
OLD-F-C-23_S12.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FYC FAC

### Female NTS, young vs aged
majiq deltapsi \
-grp1 Y-F-NTS-2_S5.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-4_S6.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-6_S7.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-7_S8.Aligned.sortedByCoord.out.majiq \
-grp2 OLD-F-NTS-14_S13.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-18_S14.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-20_S15.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-24_S16.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FYN FAN

### Male Control, young vs aged
majiq deltapsi \
-grp1 Y-M-C-25_S17.Aligned.sortedByCoord.out.majiq \
Y-M-C-27_S18.Aligned.sortedByCoord.out.majiq \
Y-M-C-31_S19.Aligned.sortedByCoord.out.majiq \
Y-M-C-33_S20.Aligned.sortedByCoord.out.majiq \
-grp2 OLD-M-C-37_S25.Aligned.sortedByCoord.out.majiq \
OLD-M-C-41_S26.Aligned.sortedByCoord.out.majiq \
OLD-M-C-43_S27.Aligned.sortedByCoord.out.majiq \
OLD-M-C-45_S28.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n MYC MAC

### Male NTS, young vs aged
majiq deltapsi \
-grp1 Y-M-NTS-26_S21.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-28_S22.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-30_S23.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-34_S24.Aligned.sortedByCoord.out.majiq \
-grp2 Aged-M-NTS-38_S1.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-42_S2.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-44_S3.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-46_S4.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n MYN MAN

### Control vs. NTS
### Female Young, control vs NTS
majiq deltapsi \
-grp1 Y-F-C-1_S1.Aligned.sortedByCoord.out.majiq \
Y-F-C-5_S2.Aligned.sortedByCoord.out.majiq \
Y-F-C--8_S3.Aligned.sortedByCoord.out.majiq \
Y-F-C-9_S4.Aligned.sortedByCoord.out.majiq \
-grp2 Y-F-NTS-2_S5.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-4_S6.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-6_S7.Aligned.sortedByCoord.out.majiq \
Y-F-NTS-7_S8.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FYC FYN

### Female Aged, control vs NTS
majiq deltapsi \
-grp1 OLD-F-C-13_S9.Aligned.sortedByCoord.out.majiq \
OLD-F-C-19_S10.Aligned.sortedByCoord.out.majiq \
OLD-F-C-21_S11.Aligned.sortedByCoord.out.majiq \
OLD-F-C-23_S12.Aligned.sortedByCoord.out.majiq \
-grp2 OLD-F-NTS-14_S13.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-18_S14.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-20_S15.Aligned.sortedByCoord.out.majiq \
OLD-F-NTS-24_S16.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n FAC FAN

### Male Young, control vs NTS
majiq deltapsi \
-grp1 Y-M-C-25_S17.Aligned.sortedByCoord.out.majiq \
Y-M-C-27_S18.Aligned.sortedByCoord.out.majiq \
Y-M-C-31_S19.Aligned.sortedByCoord.out.majiq \
Y-M-C-33_S20.Aligned.sortedByCoord.out.majiq \
-grp2 Y-M-NTS-26_S21.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-28_S22.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-30_S23.Aligned.sortedByCoord.out.majiq \
Y-M-NTS-34_S24.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n MYC MYN

### Male Aged, control vs NTS
majiq deltapsi \
-grp1 OLD-M-C-37_S25.Aligned.sortedByCoord.out.majiq \
OLD-M-C-41_S26.Aligned.sortedByCoord.out.majiq \
OLD-M-C-43_S27.Aligned.sortedByCoord.out.majiq \
OLD-M-C-45_S28.Aligned.sortedByCoord.out.majiq \
-grp2 Aged-M-NTS-38_S1.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-42_S2.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-44_S3.Aligned.sortedByCoord.out.majiq \
Aged-M-NTS-46_S4.Aligned.sortedByCoord.out.majiq \
-o ${outputDir} \
-j 20 -n MAC MAN
