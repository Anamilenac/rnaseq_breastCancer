#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5000M
#SBATCH --time=10:00:00
#SBATCH --job-name=index_bam
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_indexbam_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_indexbam_%j.e

# Go to working directory
# cd /data/users/$USER/breast_cancer

BAM_DIR=/data/users/$USER/breast_cancer/analysis/bam
SAMPLELIST=/data/users/$USER/breast_cancer/script/samplelist.tsv

BAM=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $5; exit}' $SAMPLELIST)

############################

module load UHTS/Analysis/samtools/1.10

samtools index $BAM
