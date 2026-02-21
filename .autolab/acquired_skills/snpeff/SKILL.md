# SnpEff Skill

## When to Use
Use for annotating variants with gene names, functional effects, and impact predictions.

## Standard Workflow
1. Build custom database from GFF/GBK + genome:
   - Create snpEff config entry
   - Place genome and genes files in data directory
   - Run: `snpEff build -gff3 custom_db`
2. Annotate VCF: `snpEff ann custom_db variants.vcf > annotated.vcf`
3. Parse annotations from ANN field in VCF

## Key Decisions
- For custom/assembled genomes: build a local SnpEff database using Prokka annotations
- Impact levels: HIGH, MODERATE, LOW, MODIFIER
- Parse ANN field: Allele|Annotation|Impact|Gene_Name|Gene_ID|...
- For E. coli experimental evolution, focus on MODERATE and HIGH impact variants
