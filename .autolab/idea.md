The goal is to quantify transcript expression levels from paired-end RNA-Seq reads (reads_1.fq.gz, reads_2.fq.gz) using the provided reference transcriptome (transcriptome.fa). Because the data is simulated, the quantification should exactly reproduce the underlying counts. The results represent a mapping from transcript IDs → read counts.

## Task
Perform transcript quantification on the provided paired-end RNA-Seq reads using the transcriptome reference. The output should be a .tsv file with the following columns:'transcript_id'	'count'.

## Data
Input data files are in: /private/tmp/bioagent_bench/20260221_090000/transcript-quant/data/
  reads_1.fq.gz (4.5 MB)
  reads_2.fq.gz (4.5 MB)
  transcriptome.fa (1.2 MB)

## Output
Save your final output file to: /private/tmp/bioagent_bench/20260221_090000/transcript-quant/results/
Use the exact column format specified in the task above.
