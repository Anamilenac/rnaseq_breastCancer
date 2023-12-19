
FT_COUNTS_DIR=/data/users/$USER/breast_cancer/analysis/feature_counts

# Cut columns, remove header, keep only the name of samples
cut -f1,7,8,9,10,11,12,13,14,15,16,17,18 "$FT_COUNTS_DIR/featurecounts.txt" \
| tail -n +2 \
| sed 's|/data/users/acastro/breast_cancer/analysis/bam/||g; s|_sorted.bam||g' \
> "$FT_COUNTS_DIR/clean_data_counts.txt"


