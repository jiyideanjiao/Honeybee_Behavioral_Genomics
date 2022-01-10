### honeybee genome resequencing

- Authors: Chao Tong
- Date: Jan-7-2022

- Project description:

### Qulity control
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic

fastqc {individual id}.fastq
trimmomatic PE *R1_001.fastq.gz *R2_001.fastq.gz -baseout trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
