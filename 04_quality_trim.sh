#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --job-name=quality_control_trim
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_qualitytrim_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_qualitytrim_%j.e

# Go to working directory
# cd /data/users/$USER

TRIMMING_DIR=/data/users/$USER/breast_cancer/analysis/trimming
QUALITY_TRIM_DIR=/data/users/$USER/breast_cancer/analysis/trimming/quality_trim

# create directory to save quality control for trimmed files
mkdir $QUALITY_TRIM_DIR

# module avail
module add UHTS/Quality_control/fastqc/0.11.9

# Perform quality Control
echo "Analyzing quality control"
fastqc --outdir $QUALITY_TRIM_DIR $TRIMMING_DIR/*.gz

# Generate multiqc
echo "Generating MultiQC report"
module add UHTS/Analysis/MultiQC/1.8
multiqc --outdir $QUALITY_TRIM_DIR $QUALITY_TRIM_DIR

