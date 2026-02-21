# Salmon - RNA-Seq Transcript Quantification

## When to Use
Use Salmon for fast, accurate transcript-level quantification from RNA-Seq reads. Salmon uses quasi-mapping to align reads directly to a transcriptome reference without requiring genome alignment.

## Standard Workflow

### 1. Build Index
```bash
salmon index -t transcriptome.fa -i salmon_index
```

### 2. Quantify (paired-end)
```bash
salmon quant -i salmon_index -l A \
  -1 reads_1.fq.gz -2 reads_2.fq.gz \
  -o salmon_output \
  --validateMappings
```
- `-l A` auto-detects library type
- `--validateMappings` improves accuracy

### 3. Output
Results in `salmon_output/quant.sf`:
- `Name`: transcript ID
- `NumReads`: estimated read count
- `TPM`: transcripts per million
- `Length`, `EffectiveLength`: transcript lengths

### Key Decisions
- For simulated data with exact counts, use `NumReads` column and round to integers
- Use `--validateMappings` for best accuracy
- For transcript-level quantification (not gene-level), no additional summarization needed

## Installation
```bash
conda install -c bioconda salmon
# or
mamba install -c bioconda salmon
```
