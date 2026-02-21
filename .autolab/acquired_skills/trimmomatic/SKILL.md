# Trimmomatic - Read Quality Trimming

## When to Use
Use Trimmomatic to trim adapter sequences and low-quality bases from Illumina sequencing reads.

## Standard Workflow
1. Install: `conda install -c bioconda trimmomatic`
2. Run: `trimmomatic PE <input_R1.fastq.gz> <input_R2.fastq.gz> <output_R1_paired.fastq.gz> <output_R1_unpaired.fastq.gz> <output_R2_paired.fastq.gz> <output_R2_unpaired.fastq.gz> ILLUMINACLIP:<adapters.fa>:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`

## Key Parameters
- PE: paired-end mode
- ILLUMINACLIP: adapter trimming (adapter_file:seed_mismatches:palindrome_clip_threshold:simple_clip_threshold)
- LEADING/TRAILING: remove leading/trailing low quality bases
- SLIDINGWINDOW: sliding window trimming (window_size:required_quality)
- MINLEN: minimum read length after trimming

## Key Decisions
- Use TruSeq3-PE.fa for TruSeq adapters
- Adjust quality thresholds based on FastQC results
