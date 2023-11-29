#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=16:00:00
#SBATCH --job-name=future_counts
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/log/output_featurecounts_%j.o
#SBATCH --error=/data/users/acastro/log/error_featurecounts_%j.e

BAM_DIR=/data/users/$USER/breast_cancer/analysis/bam
REFERENCE_DIR=/data/users/$USER/breast_cancer/data/reference
FT_COUNTS_DIR=/data/users/$USER/breast_cancer/analysis/feature_counts

# Go to working directory
# cd /data/users/$USER

module load UHTS/Analysis/subread/2.0.1

featureCounts -p -t exon -g gene_id -a \
    $REFERENCE_DIR/Homo_sapiens.GRCh38.110.gtf \
    -o $FT_COUNTS_DIR/featurecounts.txt \
    $BAM_DIR/*_sorted.bam 

