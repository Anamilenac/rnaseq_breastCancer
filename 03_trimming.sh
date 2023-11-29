#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --job-name=trimmomatic
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/log/output_trim_%j.o
#SBATCH --error=/data/users/acastro/log/error_trim_%j.e

READS_DIR=/data/users/$USER/breast_cancer/data/reads
TRIMMING_DIR=/data/users/$USER/breast_cancer/analysis/trimming

# Go to working directory
# cd /data/users/$USER

# Create directories to sava results of trimming
mkdir --parents $TRIMMING_DIR

# module avail
module add UHTS/Analysis/trimmomatic/0.36

echo "Trimming process"

# List of sample names
samples=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

echo "Trimming process"

for sample in "${samples[@]}"; do
    
    READS_DIR_R1=/$READS_DIR/${sample}_R1.fastq.gz
    READS_DIR_R2=/$READS_DIR/${sample}_R2.fastq.gz
    
    FORWARD_PAIR="$TRIMMING_DIR/${sample}_trimmed.paired_R1.fastq.gz"
    REVERSE_PAIR="$TRIMMING_DIR/${sample}_trimmed.paired_R2.fastq.gz"
    
    echo $sample
    
    trimmomatic PE -threads 1 -phred33 $READS_DIR_R1 $READS_DIR_R2 $FORWARD_PAIR $REVERSE_PAIR LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
    # LEADING:20 TRAILING:20: Removes leading and trailing bases with a Phred score below 20.
    # SLIDINGWINDOW:4:20: Performs a sliding window trimming, cutting when the average quality within the window falls below 20.
    # MINLEN:50: Discards reads below a length of 50 bases after trimming.
    
done

