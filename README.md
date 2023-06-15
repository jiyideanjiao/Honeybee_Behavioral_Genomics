#### Honeybee genome resequencing

- Authors: Chao Tong
- Date: Jun-1-2023

- Project description:

#### Qulity control

```
fastqc {id}.fastq
trimmomatic PE *R1_001.fastq.gz *R2_001.fastq.gz -baseout trimmed.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
```
#### Map Reads to Bee Genome

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
- 4. sort bam file
```
samtools sort {id}.bam -o {id}.sorted.bam
```
- 5. remove PCR duplication
```
samtools rmdup {id}.sorted.bam {id}_nopcr.bam
```
```
picard MarkDuplicates I={id}.bam O={id}.rmd.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
```
```
samtools index *.rmd.bam
```
- 6. calling variant
```
bcftools mpileup -f <ReferenceGenome> <All file names> | 
bcftools call --skip-variants indels --variants-only -mv -Oz -o output.vcf.gz
```
### Filter variant
```
gunzip output.vcf.gz
bcftools filter -o bee.snps.filtered.vcf -i 'QUAL>20 &&DP>5' output.vcf
```


- 1. extract dosage and genotype likelihoods
```
bcftools +dosage output.vcf > output.dosage.txt
```
```
bcftools +dosage behavior_filter1.vcf > behavior.dosage1
bcftools +dosage behavior_filter2.vcf > behavior.dosage2
bcftools +dosage behavior_filter3.vcf > behavior.dosage3
bcftools +dosage chc_filter2.vcf > chc.dosage2
bcftools +dosage chc_filter3.vcf > chc.dosage3
bcftools +dosage chc_filter4.vcf > chc.dosage4
bcftools +dosage chc_filter5.vcf > chc.dosage5
```
```
bcftools +dosage bee_filter8.vcf > bee.dosage8
```
#### Prepare Genotype and Phenotype files for running GEMMA

- 1. create relatedness matrix
```
gemma -g <Genotype> -p <Phenotype> -gk 1 -o <Relatedness>
```
```
gemma -g behavior.dosage1 -p behavior_phenotype.txt -gk 1 -o relatedness_behavior1
gemma -g behavior.dosage2 -p behavior_phenotype.txt -gk 1 -o relatedness_behavior2
gemma -g behavior.dosage3 -p behavior_phenotype.txt -gk 1 -o relatedness_behavior3
```
```
gemma -g bee.dosage8 -p bee_phenotypes.txt -gk 1 -o relatedness
```

- 2. run GEMMA
```
gemma -g <Genotype> -p <Phenotype> -n <num of column: Trait> -k <relatedness> -lm 4 -o <Trait name>
```

#### Prepare Genotype and Phenotype file for running GEMMA with large scale of variants
- 1. prepare the bed file
```
vcftools --vcf bee.vcf --plink --out output
plink --file output --make-bed --out final
plink --bfile final --freq --out MAF_check
```
- The first command will generate two output files: output.ped and output.map
- The second command will generate five output files: final.bed, final.bim, final.fam, final.log, final.nosex
- 2. prepare relatedness matrix
```
gemma -bfile final -p <phenotype> -gk -o <relatedness>
gemma -bfile final -p <Phenotype> -n <num of column: Trait> -k <relatedness> -lm 4 -o <Trait name>
```
#### Detect variant (optional)

- 1. generate bcf file
```
samtools mpileup -gf GCF_003254395.2_Amel_HAv3.1_genomic.fna {id}_nopcr.bam > {id}.bcf
```
- 2. variantion detection
```
bcftools call -vm {id}.bcf -o {id}.variants.bcf
bcftools view -v snps,indels {id}.variants.bcf > {id}.snps.vcf
```

#### Annotate Variant With Annovar
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

#### Annotate Variant With snpEff

