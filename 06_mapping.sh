#!/usr/bin/env bash

#SBATCH --array=1-12
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=13:00:00
#SBATCH --job-name=slurm_array_mapping
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_mapping_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_mapping_%j.e
#SBATCH --partition=pall

# Go to working directory
# cd /data/users/$USER/breast_cancer

REFERENCE_DIR=/data/users/$USER/breast_cancer/data/reference
BAM_DIR=/data/users/$USER/breast_cancer/analysis/bam
SAMPLELIST=/data/users/$USER/breast_cancer/script/samplelist.tsv

# Create directory to save sam files
mkdir --parents $BAM_DIR

SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

SAM_FILE="$BAM_DIR/${SAMPLE}.sam"

############################

module load UHTS/Aligner/hisat/2.2.1

# Align with Hisat 2
hisat2 -x $REFERENCE_DIR/hisat2/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
-1 $READ1 \
-2 $READ2 \
-S $SAM_FILE
