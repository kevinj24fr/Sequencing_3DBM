#!/bin/bash

#####################################################################################################
# Script to generate count matrices out of illumina generated fastq files                           #
# Author: Kevin Joseph (kevin.joseph@uniklinik-freiburg.de)                                         #
# 3DBM & NE, Neurosurgery, University Hospital of Freiburg                                          #
#####################################################################################################

echo "Setting up to process illumina sequencing data"
echo "This script is for MacOS based systems"

# Check if sufficient arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path_to_fastq_folder> <path_to_reference_genome>"
    exit 1
fi

# Assign variables
FASTQ_FOLDER=$1
REFERENCE_GENOME=$2
THREADS=4  # Set the number of threads to use

# Paths to reference files
HISAT_INDEX_BASE="${REFERENCE_GENOME}_hisat_index"

# Step 1: Create HISAT2 Index if not already present
echo "Checking HISAT2 index..."
if [ ! -f "${HISAT_INDEX_BASE}.1.ht2" ]; then
    echo "Creating HISAT2 index from reference genome..."
    hisat2-build "$REFERENCE_GENOME" "$HISAT_INDEX_BASE"
fi

# Step 2: Quality Control (FastQC)
echo "Running FastQC on raw reads..."
for FILE in "$FASTQ_FOLDER"/*.fastq.gz; do
    fastqc "$FILE" -o "$FASTQ_FOLDER"/fastqc_reports
done

# Step 3: Trimming adapters and low-quality reads with Trimmomatic
echo "Trimming reads with Trimmomatic..."
for FILE in "$FASTQ_FOLDER"/*_R1_001.fastq.gz; do
    BASE=$(basename "$FILE" "_R1_001.fastq.gz")
    R1="$FASTQ_FOLDER/${BASE}_R1_001.fastq.gz"
    R2="$FASTQ_FOLDER/${BASE}_R2_001.fastq.gz"
    TRIMMED_R1="$FASTQ_FOLDER/${BASE}_R1_trimmed.fastq"
    TRIMMED_R2="$FASTQ_FOLDER/${BASE}_R2_trimmed.fastq"
    
    trimmomatic PE -threads $THREADS "$R1" "$R2" "$TRIMMED_R1" "${TRIMMED_R1}_unpaired" \
    "$TRIMMED_R2" "${TRIMMED_R2}_unpaired" SLIDINGWINDOW:4:20 MINLEN:36
done

# Step 4: Aligning with HISAT2
echo "Aligning reads with HISAT2..."
for FILE in "$FASTQ_FOLDER"/*_R1_trimmed.fastq; do
    BASE=$(basename "$FILE" "_R1_trimmed.fastq")
    TRIMMED_R1="$FASTQ_FOLDER/${BASE}_R1_trimmed.fastq"
    TRIMMED_R2="$FASTQ_FOLDER/${BASE}_R2_trimmed.fastq"
    SAM_FILE="$FASTQ_FOLDER/${BASE}.sam"
    
    hisat2 -p $THREADS -x "$HISAT_INDEX_BASE" -1 "$TRIMMED_R1" -2 "$TRIMMED_R2" -S "$SAM_FILE"
done

# Step 5: Convert SAM to BAM and sort using Samtools
echo "Converting SAM to BAM and sorting..."
for FILE in "$FASTQ_FOLDER"/*.sam; do
    BASE=$(basename "$FILE" ".sam")
    BAM_FILE="$FASTQ_FOLDER/${BASE}.bam"
    SORTED_BAM_FILE="$FASTQ_FOLDER/${BASE}_sorted.bam"
    
    samtools view -bS "$FILE" > "$BAM_FILE"
    samtools sort "$BAM_FILE" -o "$SORTED_BAM_FILE"
    samtools index "$SORTED_BAM_FILE"
done

# Step 6: Count reads with featureCounts
echo "Counting reads with featureCounts..."
BAM_FILES=$(ls "$FASTQ_FOLDER"/*_sorted.bam | tr '\n' ' ')
featureCounts -T $THREADS -a "$REFERENCE_GENOME".gtf -o "$FASTQ_FOLDER"/counts.txt $BAM_FILES

echo "Pipeline complete. Count matrix saved to $FASTQ_FOLDER/counts.txt"
