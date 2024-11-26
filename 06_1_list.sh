#!/bin/bash

READS_DIR=/data/users/$USER/breast_cancer/data/reads
SAMPLE_LIST=/data/users/$USER/breast_cancer/script/
BAM_DIR=/data/users/$USER/breast_cancer/analysis/bam

if [ -f samplelist.tsv ]; then
    echo "Remove file"
    rm samplelist.tsv
fi

for FILE in $READS_DIR/*_*1.fastq.gz; do
    PREFIX="${FILE%_*.fastq.gz}"
    SAMPLE=$(basename $PREFIX)
    SAM_FILE="${BAM_DIR}/${SAMPLE}.sam"
    BAM_FILE="${BAM_DIR}/${SAMPLE}_sorted.bam"
    echo -e "${SAMPLE}\t$FILE\t${FILE%?.fastq.gz}2.fastq.gz\t${SAM_FILE}\t${BAM_FILE}" >>samplelist.tsv
done
