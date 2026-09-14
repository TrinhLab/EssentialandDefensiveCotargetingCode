#!/usr/bin/env bash
# Run three short-sequence BLAST searches against a pre-built database.
#
# Inputs:
#   1. BLAST database basename
#   2. Output directory
#   3-5. Three guide FASTA files, one for each batch
#
# Outputs:
#   One tabular BLAST result file per batch. The columns are query ID, subject
#   ID, percent identity, alignment length, mismatches, gap openings, query and
#   subject coordinates, E-value, and bit score.

set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 DB_NAME OUTPUT_DIR GUIDE_BATCH_1 GUIDE_BATCH_2 GUIDE_BATCH_3" >&2
    exit 2
fi

DB_NAME=$1
OUTPUT_DIR=$2
GUIDE_FILES=("${@:3}")

mkdir -p "$OUTPUT_DIR"

echo "Starting 3 guide-search batches."
echo "BLAST database: $DB_NAME"
echo "Output directory: $OUTPUT_DIR"

for batch_index in "${!GUIDE_FILES[@]}"; do
    batch=$((batch_index + 1))
    guide_file=${GUIDE_FILES[$batch_index]}
    output_file="$OUTPUT_DIR/guide_blast_results_${batch}.tsv"

    if [[ ! -f "$guide_file" ]]; then
        echo "Missing guide input: $guide_file" >&2
        exit 1
    fi

    echo "Running batch $batch of 3: $guide_file"
    blastn \
        -task blastn-short \
        -query "$guide_file" \
        -db "$DB_NAME" \
        -out "$output_file" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -num_threads "${BLAST_THREADS:-1}" \
        -perc_identity 65 \
        -qcov_hsp_perc 65 \
        -max_target_seqs 100000000 \
        -evalue 10000000

    echo "Batch $batch complete: $output_file"
done

echo "All 3 batches completed."
