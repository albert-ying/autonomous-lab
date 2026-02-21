# Metagenomics Analysis Pipeline

## When to Use
Use for end-to-end metagenomic community profiling from raw reads to taxonomic abundance tables.

## Standard Workflow
1. Quality control: Trim adapters and low-quality bases with Trimmomatic
2. Taxonomic classification: Classify reads with Kraken2 against reference database
3. Abundance estimation: Estimate taxon abundances with Bracken
4. Parse results: Generate summary tables (CSV) with taxonomic assignments and relative abundances

## Pipeline Steps
1. QC: trimmomatic PE -> paired reads
2. Classification: kraken2 --db <db> --paired -> kraken report
3. Abundance: bracken -d <db> -i kraken_report -l P -> phylum-level abundances
4. Parse: extract OTU, Kingdom, Phylum, and relative abundances per sample

## Key Decisions
- Choose appropriate taxonomic level for Bracken (P for phylum)
- Determine read length from data to select correct kmer_distrib file
- Merge multi-sample results into single table with relative abundances
