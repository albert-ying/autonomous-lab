# Bracken - Bayesian Reestimation of Abundance with KrakEN

## When to Use
Use Bracken to estimate species/genus/phylum-level abundances from Kraken2 classification results.

## Standard Workflow
1. Install: `conda install -c bioconda bracken`
2. Run after Kraken2: `bracken -d <db_path> -i <kraken_report.txt> -o <bracken_output.txt> -r <read_length> -l <level> -t <threshold>`
3. Output: abundance table with fraction_total_reads per taxon

## Key Parameters
- `-d`: Kraken2 database path (must contain kmer_distrib files)
- `-i`: Kraken2 report file
- `-o`: Output file
- `-r`: Read length (50, 75, 100, 150, 200, 250, 300)
- `-l`: Taxonomic level (S=species, G=genus, F=family, O=order, C=class, P=phylum, K=kingdom, D=domain)
- `-t`: Minimum read threshold (default 0)

## Key Decisions
- Read length should match the kmer_distrib files available in the database
- Taxonomic level determines granularity of abundance estimates
- The kmer_distrib files (e.g., database150mers.kmer_distrib) are required
