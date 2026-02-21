# Prodigal Skill

## When to Use
Use Prodigal for ab initio gene prediction in prokaryotic genomes.

## Standard Workflow
1. Run Prodigal on each genome FASTA: `prodigal -i genome.fna -a proteins.faa -o genes.gff -f gff`
2. Extract protein sequences from -a output
3. Extract CDS coordinates from GFF output

## Key Decisions
- Use `-p single` for single genome mode (default)
- Use `-p meta` for metagenomics mode
- Output both protein (-a) and nucleotide (-d) sequences
