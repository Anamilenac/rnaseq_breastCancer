#!/usr/bin/env bash

PATH_DATA=/data/users/acastro/rnaseq_course/breast_cancer/data/
SOURCE=/data/courses/rnaseq_course/breastcancer_de

# rm -r $PATH_DATA/reads

echo "copying reads"
cp --verbose --recursive $SOURCE $PATH_DATA

# echo "uncompressing reads"
# gunzip -k $PATH_DATA/breastcancer_de/reads/*.gz


