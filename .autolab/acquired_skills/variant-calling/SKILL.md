# Variant Calling Skill

## When to Use
Use when calling SNPs and indels from aligned BAM files against a reference.

## Standard Workflow
1. Mark duplicates (optional): `samtools markdup`
2. Call variants with freebayes: `freebayes -f reference.fasta -p 1 sample.bam > variants.vcf`
   OR with bcftools: `bcftools mpileup -f ref.fa sample.bam | bcftools call -mv -Oz -o variants.vcf.gz`
3. Filter variants: `bcftools filter -s LowQual -e 'QUAL<20' variants.vcf`

## Key Decisions
- For haploid organisms (E. coli), use ploidy=1
- freebayes is good for bacterial genomes; bcftools is fast and reliable
- Filter on QUAL score and read depth
