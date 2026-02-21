# OrthoFinder Skill

## When to Use
Use OrthoFinder for identifying orthologous gene groups across multiple genomes, inferring gene trees, and identifying orthogroups (clusters of orthologous genes).

## Standard Workflow
1. Prepare protein FASTA files (one per genome) from gene predictions
2. Run OrthoFinder: `orthofinder -f <protein_fasta_dir>`
3. Parse results from Orthogroups/Orthogroups.tsv
4. Filter orthogroups present in all species (core orthogroups)

## Key Decisions
- Use protein sequences (not nucleotide) for ortholog detection
- Choose appropriate inflation parameter for MCL clustering
- Filter single-copy orthologs for phylogenetic analysis
