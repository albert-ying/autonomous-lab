# SPAdes Assembly Skill

## When to Use
Use for de novo genome assembly when no reference genome is available.

## Standard Workflow
1. Run SPAdes: `spades.py -1 R1.fastq.gz -2 R2.fastq.gz -o assembly_output --careful`
2. Check assembly stats: look at scaffolds.fasta or contigs.fasta
3. Use assembled genome as reference for read mapping

## Key Decisions
- Use `--careful` flag for bacterial genomes to reduce misassemblies
- For small bacterial genomes, default k-mer sizes work well
- Output scaffolds.fasta is typically preferred over contigs.fasta
- Assembly quality can be checked with QUAST
- For experimental evolution: assemble ancestor, use as reference for evolved lines
