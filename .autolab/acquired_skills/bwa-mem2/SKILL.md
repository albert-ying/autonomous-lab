# BWA-MEM2 Skill

## When to Use
Use for aligning paired-end Illumina reads to a reference genome. BWA-MEM2 is the successor to BWA-MEM with improved performance.

## Standard Workflow
1. Index reference: `bwa index reference.fasta`
2. Align reads: `bwa mem -t 4 reference.fasta R1.fastq.gz R2.fastq.gz | samtools sort -o aligned.bam`
3. Index BAM: `samtools index aligned.bam`

## Key Decisions
- Use `-t` flag to set number of threads
- Pipe directly to samtools sort for efficiency
- Add read group info with `-R` flag: `bwa mem -R '@RG\tID:sample\tSM:sample' ...`
- For bacterial genomes, default parameters work well
