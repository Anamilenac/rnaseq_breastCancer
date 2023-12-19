#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --job-name=indexing
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/breast_cancer/log/output_indexing_%j.o
#SBATCH --error=/data/users/acastro/breast_cancer/log/error_indexing_%j.e

# Go to working directory
# cd /data/users/$USER

REFERENCE_DIR=/data/users/$USER/breast_cancer/data/reference
HISAT2_DIR=/data/users/$USER/breast_cancer/data/reference/hisat2

# Create directory to save reference genoma
mkdir --parents $REFERENCE_DIR

# Create directory to save index
mkdir --parents $HISAT2_DIR

# Copy locally reference genome and index
cp --verbose --recursive /data/courses/rnaseq_course/breastcancer_de/shared_workspace/Homo_sapiens.GRCh38.110.gtf $REFERENCE_DIR
cp --verbose --recursive /data/courses/rnaseq_course/breastcancer_de/Homo_sapiens.GRCh38.dna.primary_assembly.fa $REFERENCE_DIR

# Optional, download the reference genome
# wget -P $REFERENCE_DIR https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# module avail
module add UHTS/Aligner/hisat/2.2.1

# # Unconpress reads
# echo "uncompressing reference"
# gunzip $REFERENCE_DIR/*.fa.gz

echo "Indexing the reference genoma"
hisat2-build -p 16 $REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa $HISAT2_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa
















