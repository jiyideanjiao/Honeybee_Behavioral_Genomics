#### Honeybee genome resequencing

- Authors: Chao Tong
- Date: Jun-1-2023

- Project description:

#### Qulity control
```
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
```
```
fastqc {id}.fastq
trimmomatic PE *R1_001.fastq.gz *R2_001.fastq.gz -baseout trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
```
#### Mapping Reads to Bee Genome
```
conda install -c bioconda bwa
```
- 1. build genome index
```
bwa index GCF_003254395.2_Amel_HAv3.1_genomic.fna
```
- 2. map reads to honeybee genome to generate sam file
```
bwa mem GCF_003254395.2_Amel_HAv3.1_genomic.fna {id}_R1.fq {id}_R2.fq > {id}.sam
```
- 3. transform sam file to bam file
```
samtools faidx GCF_003254395.2_Amel_HAv3.1_genomic.fna
samtools view -bhS -t GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai {id}.sam -o {id}.bam
```
- 4. sorting bam file
```
samtools sort {id}.bam -o {id}.sorted.bam
```
- 5. build index for bam file
```
samtools index {id}.sorted.bam
```
- 6. estimate the depth
```
samtools depth {id}.sorted.bam > {id}_depth.txt
samtools flagstat {id}.sorted.bam > {id}_counts.txt
```

#### Preparing Genotype and Phenotype file for running GEMMA

- 1. remove PCR duplication
```
samtools rmdup {id}.sorted.bam {id}_nopcr.bam
```
- 2. variant calling
option 1
```
samtools mpileup -g -f <ReferenceGenome> <All file names> | 
bcftools call --skip-variants indels --variants-only -mv -Oz > output.vcf.gz
```
option 2
```
bcftools mpileup -f <ReferenceGenome> <All file names> | 
bcftools call --skip-variants indels --variants-only -mv -Oz > output.vcf.gz
```
- 3. variant filtering
```
gunzip output.vcf.gz
bcftools filter -o bee.snps.filtered.vcf -i 'QUAL>20 &&DP>5' output.vcf
```
- 4. extract dosage and genotype likelihoods
```
bcftools +dosage output.vcf > output.dosage.txt
```
- 5. create relatedness matrix
option 1
```
gemma -g Genotype.txt -p Phenotype.txt -gk 1 -o Relatedness.txt
```
option 2
```
vcftools --vcf bee.snps.filtered.vcf --plink --out output
```
This command will generate two output files: output.ped and output.map
```
vcftools --ped output.ped --map output.map --relatedness2 --out relatedness
```

- 5. GEMMA run
```
gemma -g Genotype.txt -p Phenotype.txt -n <num of trait> -k Relatedness.txt -lm 4 -o Results.txt
```

#### Detect variant (optional)
- 1. remove PCR duplication
```
samtools rmdup {id}.sorted.bam {id}_nopcr.bam
```
- 2. generate bcf file
```
samtools mpileup -gf GCF_003254395.2_Amel_HAv3.1_genomic.fna {id}_nopcr.bam > {id}.bcf
```
- 3. variantion detection
```
bcftools call -vm {id}.bcf -o {id}.variants.bcf
bcftools view -v snps,indels {id}.variants.bcf > {id}.snps.vcf
```




#### Annotate variant
- 1. generate annovar input file
download and install **Annovar** [link](http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz)
```
convert2annovar.pl -format vcf4 {id}.snps.vcf > {id}.snps.avinput
```
- 2. build local bee database
```
conda install -c bioconda ucsc-gff3togenepred
```
```
gff3ToGenePred -useName GCF_003254395.2_Amel_HAv3.1_genomic.gff genome_refGene.txt
```
```
cut -f 12 genome_refGene.txt > column1.txt
cut -f 2-15 genome_refGene.txt > column_else.txt
paste column1.txt column_else.txt > genome_new_refGene.txt
```
```
retrieve_seq_from_fasta.pl -format refGene -seqfile GCF_003254395.2_Amel_HAv3.1_genomic.fna -outfile genome_new_refGeneMrna.fa genome_new_refGene.txt
cp genome_new_refGene* ./beedb/
```
```
annotate_variation.pl --geneanno --dbtype refGene --buildver genome_new  {id}.snps.avinput ./beedb/
mv *_function result
```
