#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --job-name=run_trimmomatic
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/rnaseq_course/output_fastqc_%j.o
#SBATCH --error=/data/users/acastro/rnaseq_course/error_fastqc_%j.e

# module avail
module add UHTS/Analysis/trimmomatic/0.36

# rm --recursive $PATH_OUTPUT_TRIM

# List of sample names
samples=("HER21" "HER22" "HER23" "NonTNBC1" "NonTNBC2" "NonTNBC3" "Normal1" "Normal2" "Normal3" "TNBC1" "TNBC2" "TNBC3")

PATH_OUTPUT_TRIM="/data/users/acastro/rnaseq_course/breast_cancer/analysis_step/trimmed_reads"

mkdir --parents "$PATH_OUTPUT_TRIM"

echo "Trimming process"


for sample in "${samples[@]}"; do

    PATH_READ1="/data/users/acastro/rnaseq_course/breast_cancer/data/breastcancer_de/reads/${sample}_R1.fastq.gz"
    PATH_READ2="/data/users/acastro/rnaseq_course/breast_cancer/data/breastcancer_de/reads/${sample}_R2.fastq.gz"

    OUTPUT_FORWARD_PAIR="$PATH_OUTPUT_TRIM/${sample}_trimmed.paired_R1.fastq.gz"
    OUTPUT_FORWARD_UNPAIR="$PATH_OUTPUT_TRIM/${sample}_trimmed.unpaired_R1.fastq.gz"
    OUTPUT_REVERSE_PAIR="$PATH_OUTPUT_TRIM/${sample}_trimmed.paired_R2.fastq.gz"
    OUTPUT_REVERSE_UNPAIR="$PATH_OUTPUT_TRIM/${sample}_trimmed.unpaired_R2.fastq.gz"

    echo $sample

    trimmomatic PE -threads 1 -phred33 \
    $PATH_READ1 $PATH_READ2 \
    $OUTPUT_FORWARD_PAIR $OUTPUT_FORWARD_UNPAIR \
    $OUTPUT_REVERSE_PAIR $OUTPUT_REVERSE_UNPAIR \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

done


