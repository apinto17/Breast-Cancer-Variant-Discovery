#!/bin/bash

# WARNING! COMPUTATIONALLY INTENSIVE!! 
# The resulting fastq files are 17GB, must of this script had to be run overnight, the whole thing took multiple days to run
# If you're interested in running this, try running it on a cluster

# getting data
fasterq-dump --split-files --threads 4 SRR925782.sra
fastqc SRR925782_1.fastq

# trimming data
# Did not pass base read quality, dropped off at the end
fastp -i SRR925782_1.fastq -o SRR925782_1.trimmed.fastq -h fastp.html -j fastp.json 

# left over some adapter sequences
fastp -i SRR925782_1.fastq -o SRR925782_1.cleaned.fastq --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 -h fastp.html 

# trimmed adapters
cutadapt \
  -a TGGAATTCTCGGGTGCCAAGG \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -q 20 \
  -m 50 \
  -o SRR925782_1.cleaned.v2.fastq \
  SRR925782_1.cleaned.fastq

# WARNING! Computationally intensive!
# index reference genome
bwa index GRCh38.p14.genome.fa

# WARNING! Computationally intensive!
# aligned to reference
bwa mem -t 8 GRCh38.p14.genome.fa SRR925782_1.cleaned.v2.fastq > SRR925782_1.aligned.sam

# clean up sam file, last line was corrupted
sed '$d' SRR925782_1.aligned.sam > SRR925782_1.aligned.clean.sam

# WARNING! Computationally intensive!
# convert to bam file
samtools view -F 4 -b SRR925782_1.aligned.clean.sam > aligned.bam

# sort
samtools sort aligned.bam -o aligned_sorted.bam

# index
samtools index aligned_sorted.bam


# mark duplicates
# need to install picard.jar
java -jar picard.jar MarkDuplicates \
  I=aligned_sorted.bam \
  O=aligned_sorted.markdup.bam \
  M=marked_dup_metrics.txt \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true

# WARNING! Computationally intensive!
# index reference
samtools faidx GRCh38.p14.genome.fa


# tried gatk, but got read groups error, so used Picard to add read groups
picard AddOrReplaceReadGroups \
  I= aligned_sorted.markdup.bam \
  O=rg_added.bam \
  RGID=1 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=sample1


# re-sorted and indexed
samtools sort -o rg_added.sorted.bam rg_added.bam
samtools index rg_added.sorted.bam

# then ran Gatk to find variants
gatk HaplotypeCaller \
  -R GRCh38.p14.genome.fa \
  -I rg_added.sorted.bam \
  -O output.vcf.gz