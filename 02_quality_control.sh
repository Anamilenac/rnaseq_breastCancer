#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --job-name=quality_control
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_quality_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_quality_%j.e

# Go to working directory
# cd /data/users/$USER

READS_DIR=/data/users/$USER/breast_cancer/data/reads
QUALITY_CONTROL_DIR=/data/users/$USER/breast_cancer/analysis/quality_control

# Create directories to results of quality control
mkdir --parents $QUALITY_CONTROL_DIR

# module avail
module add UHTS/Quality_control/fastqc/0.11.9

# Perform quality Control
echo "Analyzing quality control"
fastqc --outdir $QUALITY_CONTROL_DIR $READS_DIR/*.gz

# Generate multiqc
echo "Generating MultiQC report"
module add UHTS/Analysis/MultiQC/1.8
multiqc --outdir $QUALITY_CONTROL_DIR $QUALITY_CONTROL_DIR


