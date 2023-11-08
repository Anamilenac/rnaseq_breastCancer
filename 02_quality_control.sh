#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --job-name=run_fastqc
#SBATCH --mail-user=ana.castromarquez@students.unibe.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/data/users/acastro/slurm_tutorial/output_fastqc_%j.o
#SBATCH --error=/data/users/acastro/slurm_tutorial/error_fastqc_%j.e

# module avail
module add UHTS/Quality_control/fastqc/0.11.9

PATH_READS=/data/users/acastro/rnaseq_course/breast_cancer/data/breastcancer_de/reads
PATH_OUTPUT_ANLYSIS=/data/users/acastro/rnaseq_course/breast_cancer/analysis_step/quality_control

rm --recursive $PATH_OUTPUT_ANLYSIS

mkdir --parents $PATH_OUTPUT_ANLYSIS 

echo "Analyzing quality control"
fastqc --outdir $PATH_OUTPUT_ANLYSIS $PATH_READS/*.gz
