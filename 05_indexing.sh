#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --job-name=indexing
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/log/output_indexing_%j.o
#SBATCH --error=/data/users/acastro/log/error_indexing_%j.e

REFERENCE_DIR=/data/users/$USER/breast_cancer/data/reference
HISAT2_DIR=/data/users/$USER/breast_cancer/data/reference/hisat2

# Go to working directory
# cd /data/users/$USER

Create directory to save reference genoma 
# mkdir --parents $REFERENCE_DIR

# Create directory to save index
mkdir --parents $HISAT2_DIR

# Download the reference genome
# wget -P $REFERENCE_DIR https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Indexing the reference genoma

# module avail
module add UHTS/Aligner/hisat/2.2.1

# # Unconpress reads
# echo "uncompressing reference"
# gunzip $REFERENCE_DIR/*.fa.gz

echo "Indexing the reference genoma"
hisat2-build -p 16 $REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa $HISAT2_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# optional instead of download the reference genome
# cp --verbose --recursive /data/courses/rnaseq_course/breastcancer_de/shared_wor $READS_DIR













