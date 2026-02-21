# BWA Mapping Skill

## When to Use
Use when mapping paired-end or single-end sequencing reads to a reference genome.

## Standard Workflow
1. Index reference genome: `bwa index reference.fasta`
2. Map reads: `bwa mem -t <threads> reference.fasta R1.fastq.gz R2.fastq.gz > aligned.sam`
3. Convert to sorted BAM: `samtools sort -o aligned.sorted.bam aligned.sam`
4. Index BAM: `samtools index aligned.sorted.bam`

## Key Decisions
- Use BWA-MEM for reads >70bp (standard for Illumina HiSeq)
- Add read group info with -R flag for downstream variant calling
- Use appropriate thread count for available CPU cores
