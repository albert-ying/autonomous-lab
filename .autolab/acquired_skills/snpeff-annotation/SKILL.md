# SnpEff Annotation Skill

## When to Use
Use when annotating variants with functional impact predictions (gene, effect, impact).

## Standard Workflow
1. Build SnpEff database from reference genome + GFF: configure snpEff.config, then `snpEff build -gff3 -v genome_name`
2. Annotate VCF: `snpEff ann genome_name variants.vcf > annotated.vcf`
3. Parse annotations from ANN field in VCF

## Key Decisions
- For de novo assembled genomes, build a custom SnpEff database using Prokka annotations
- Impact levels: HIGH, MODERATE, LOW, MODIFIER
- Parse ANN field: Allele|Annotation|Impact|Gene_Name|...
