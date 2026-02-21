# Kallisto - RNA-Seq Pseudoalignment Quantification

## When to Use
Use Kallisto as an alternative to Salmon for transcript-level quantification from RNA-Seq reads using pseudoalignment.

## Standard Workflow

### 1. Build Index
```bash
kallisto index -i kallisto_index transcriptome.fa
```

### 2. Quantify (paired-end)
```bash
kallisto quant -i kallisto_index -o kallisto_output reads_1.fq.gz reads_2.fq.gz
```

### 3. Output
Results in `kallisto_output/abundance.tsv`:
- `target_id`: transcript ID
- `est_counts`: estimated read count
- `tpm`: transcripts per million

## Installation
```bash
conda install -c bioconda kallisto
```
