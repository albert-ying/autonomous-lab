# Comparative Genomics Skill

## When to Use
Use for comparing gene content, synteny, and functional modules across genomes.

## Standard Workflow
1. Predict genes with Prodigal
2. Identify orthologous clusters with OrthoFinder or similar
3. Annotate clusters with functional databases (COG, KEGG, etc.)
4. Filter for conserved gene clusters present across all genomes
5. Identify co-evolving gene modules

## Key Decisions
- Use KEGG/COG annotations for functional characterization
- Filter for clusters present in all genomes for core analysis
- Use GFF annotations when available for annotation transfer
