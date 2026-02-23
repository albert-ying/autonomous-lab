Identify differentialy expressed genes between planktonic and biofilm conditions of Candida parapsilosis. The output should be a CSV file with the following columns: gene_id,log2FoldChange,pvalue,padj
CPAR2_00000,1.0, 2.0, 1e-44

## Data
Input data files are in: ./data/
  SRR1278968_1.fastq (3386.3 MB)
  SRR1278968_2.fastq (3386.3 MB)
  SRR1278969_1.fastq (3389.1 MB)
  SRR1278969_2.fastq (3389.1 MB)
  SRR1278970_1.fastq (3319.7 MB)
  SRR1278970_2.fastq (3319.7 MB)
  SRR1278971_1.fastq (3243.4 MB)
  SRR1278971_2.fastq (3243.4 MB)
  SRR1278972_1.fastq (3374.7 MB)
  SRR1278972_2.fastq (3374.7 MB)
  SRR1278973_1.fastq (3296.7 MB)
  SRR1278973_2.fastq (3296.7 MB)

Reference files (genome, annotations, databases) are in: ./reference/
  C_parapsilosis_CDC317_current_chromosomes.fasta (12.6 MB)
  C_parapsilosis_CDC317_current_features.gff (5.6 MB)

## Background
The dataset consists of RNA-Seq samples from Candida parapsilosis wild-type (WT) strains grown in planktonic and biofilm conditions, generated as part of a study on gene expression and biofilm formation. The samples were sequenced on the Illumina HiSeq 2000 platform. The goal of this analysis is to perform differential expression analysis using DESeq2 to identify genes that are significantly up- or down-regulated between planktonic and biofilm conditions, providing insights into biofilm-associated transcriptional changes.
