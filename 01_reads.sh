#!/usr/bin/env bash

# # Go to working directory 
# cd /data/users/$USER

# # Create directories to save reads 
# mkdir --parents /data/users/$USER/breast_cancer/data

# READS_DIR=/data/users/$USER/breast_cancer/data/reads

# #Copy reads locally 
# echo "copying reads"
# cp --verbose --recursive /data/courses/rnaseq_course/breastcancer_de/reads/ $READS_DIR

# Unconpress reads
echo "uncompressing reads"
gunzip -k $READS_DIR/*.gz


