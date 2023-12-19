#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=13:00:00
#SBATCH --job-name=ConvertSamToBam
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_samtobam_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_samtobam_%j.e

# Go to working directory
# cd /data/users/$USER

BAM_DIR=/data/users/$USER/breast_cancer/analysis/bam
SAMPLELIST=/data/users/$USER/breast_cancer/script/samplelist.tsv

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
SAM=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $4; exit}' $SAMPLELIST`

BAM_FILE="$BAM_DIR/${SAMPLE}_sorted.bam"

############################

module load UHTS/Analysis/samtools/1.10

samtools view -hbS $SAM | samtools sort -o $BAM_FILE


