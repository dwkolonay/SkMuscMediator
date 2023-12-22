#!/bin/bash



# Define input and output directories

input_dir="path/to/raw_data"

output_dir="path/to/analysis_results"



# Create output directory if it doesn't exist

mkdir -p "$output_dir"



# FastQC - Quality assessment

fastqc -o "$output_dir/fastqc_reports" "$input_dir/*.fastq"



# HISAT2 - Alignment

hisat2 -x path/to/hisat2_index -U "$input_dir/*.fastq" -S "$output_dir/aligned_reads.sam"



# Convert SAM to BAM

samtools view -bS "$output_dir/aligned_reads.sam" -o "$output_dir/aligned_reads.bam"



# Sort BAM file

samtools sort "$output_dir/aligned_reads.bam" -o "$output_dir/aligned_reads_sorted.bam"



# Index BAM file

samtools index "$output_dir/aligned_reads_sorted.bam"



# featureCounts - Read counts

featureCounts -a path/to/annotation.gtf -o "$output_dir/read_counts.txt" "$output_dir/aligned_reads_sorted.bam"



# Clean up intermediate files if needed

# rm "$output_dir/aligned_reads.sam" "$output_dir/aligned_reads.bam"



echo "Pipeline completed successfully."


