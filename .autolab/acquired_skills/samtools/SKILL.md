# Samtools Skill

## When to Use
Use for BAM/SAM file manipulation, sorting, indexing, and basic variant calling support.

## Standard Workflow
1. Convert SAM to BAM: `samtools view -bS input.sam -o output.bam`
2. Sort BAM: `samtools sort -o sorted.bam input.bam`
3. Index BAM: `samtools index sorted.bam`
4. View stats: `samtools flagstat sorted.bam`
5. Index FASTA: `samtools faidx reference.fasta`

## Key Decisions
- Always sort and index BAMs before variant calling
- Use samtools flagstat to check mapping quality
- samtools depth for coverage analysis
