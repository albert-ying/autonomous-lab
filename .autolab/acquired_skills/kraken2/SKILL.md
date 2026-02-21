# Kraken2 - Taxonomic Classification

## When to Use
Use Kraken2 for taxonomic classification of metagenomic sequencing reads against a reference database.

## Standard Workflow
1. Install: `conda install -c bioconda kraken2` or `brew install kraken2`
2. Run classification: `kraken2 --db <db_path> --paired <R1.fastq.gz> <R2.fastq.gz> --output <output.txt> --report <report.txt> --threads <N>`
3. Key output files:
   - output.txt: per-read classification
   - report.txt: summary report with read counts per taxon

## Key Parameters
- `--db`: path to Kraken2 database (directory containing hash.k2d, opts.k2d, taxo.k2d)
- `--paired`: for paired-end reads
- `--report`: generate summary report
- `--confidence`: confidence threshold (0-1, default 0)
- `--threads`: number of threads

## Key Decisions
- Confidence threshold affects sensitivity vs specificity tradeoff
- Database choice determines which taxa can be identified
