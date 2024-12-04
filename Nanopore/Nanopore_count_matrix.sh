#!/bin/bash

#####################################################################################################
# Script to generate count matrices out of illumina generated fastq files, Nanopore Edition         #
# Author: Kevin Joseph (kevin.joseph@uniklinik-freiburg.de)                                         #
# 3DBM & NE, Neurosurgery, University Hospital of Freiburg                                          #
#####################################################################################################

# Check if sufficient arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path_to_fastq_folder> <path_to_reference_genome>"
    exit 1
fi

# Assign variables
FASTQ_FOLDER="$1"
REFERENCE_GENOME="$2"
THREADS=4  # Set the number of threads to use

# Paths to reference files
MINIMAP_INDEX_BASE="${REFERENCE_GENOME}.mmi"

# Step 1: Create Minimap2 Index if not already present
echo "Checking Minimap2 index..."
if [ ! -f "$MINIMAP_INDEX_BASE" ]; then
    echo "Creating Minimap2 index from reference genome..."
    minimap2 -d "$MINIMAP_INDEX_BASE" "$REFERENCE_GENOME"
fi

# Step 2: Quality Control (FastQC)
echo "Running FastQC on raw reads..."
FASTQC_OUTPUT_DIR="$FASTQ_FOLDER/fastqc_reports"
mkdir -p "$FASTQC_OUTPUT_DIR"
for FILE in "$FASTQ_FOLDER"/*.fastq; do
    if [ -f "$FILE" ]; then
        fastqc "$FILE" -o "$FASTQC_OUTPUT_DIR"
    else
        echo "Warning: FASTQ file $FILE not found or not readable. Skipping."
    fi
done

# Step 3: Basecalling (if needed, using Guppy)
echo "Basecalling using Guppy (if applicable)..."
# Note: Adjust command based on your requirement. Skip if already basecalled.
for FILE in "$FASTQ_FOLDER"/*.fast5; do
    if [ -f "$FILE" ]; then
        echo "Basecalling $FILE with Guppy..."
        guppy_basecaller -i "$FASTQ_FOLDER" -s "$FASTQ_FOLDER/basecalled" -c dna_r9.4.1_450bps_fast.cfg -x auto --num_callers $THREADS
    fi
done

# Step 4: Aligning with Minimap2
echo "Aligning reads with Minimap2..."
ALIGN_OUTPUT_DIR="$FASTQ_FOLDER/alignment"
mkdir -p "$ALIGN_OUTPUT_DIR"
for FILE in "$FASTQ_FOLDER"/*.fastq; do
    if [ -f "$FILE" ]; then
        BASE=$(basename "$FILE" ".fastq")
        SAM_FILE="$ALIGN_OUTPUT_DIR/${BASE}.sam"
        
        minimap2 -t $THREADS -ax map-ont "$MINIMAP_INDEX_BASE" "$FILE" > "$SAM_FILE"
    else
        echo "Warning: FASTQ file $FILE not found. Skipping."
    fi
done

# Step 5: Convert SAM to BAM and sort using Samtools
echo "Converting SAM to BAM and sorting..."
for SAM_FILE in "$ALIGN_OUTPUT_DIR"/*.sam; do
    if [ -f "$SAM_FILE" ]; then
        BASE=$(basename "$SAM_FILE" ".sam")
        BAM_FILE="$ALIGN_OUTPUT_DIR/${BASE}.bam"
        SORTED_BAM_FILE="$ALIGN_OUTPUT_DIR/${BASE}_sorted.bam"
        
        samtools view -bS "$SAM_FILE" > "$BAM_FILE"
        samtools sort "$BAM_FILE" -o "$SORTED_BAM_FILE"
        samtools index "$SORTED_BAM_FILE"
    else
        echo "Warning: SAM file $SAM_FILE not found. Skipping."
    fi
done

# Step 6: Count reads with featureCounts
echo "Counting reads with featureCounts..."
BAM_FILES=$(ls "$ALIGN_OUTPUT_DIR"/*_sorted.bam 2>/dev/null | tr '\n' ' ')
if [ -z "$BAM_FILES" ]; then
    echo "Error: No BAM files found for featureCounts. Exiting."
    exit 1
fi
featureCounts -T $THREADS -a "$REFERENCE_GENOME".gtf -o "$ALIGN_OUTPUT_DIR/counts.txt" $BAM_FILES

echo "Pipeline complete. Count matrix saved to $ALIGN_OUTPUT_DIR/counts.txt"
