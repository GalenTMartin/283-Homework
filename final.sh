#!/usr/bin/env bash
#
#$ -N bismark_283
#$ -q epyc,free64,abio,bio
#$ -m beas
#$ -M galentm@uci.edu
### -o
### -e
#$ -pe openmp 10
#$ -ckpt restart

module load sratoolkit/2.8.1
#module load anaconda/3.6-5.0.1
module load bismark/0.2.1
module load trimmomatic/0.35
module load fastqc/0.11.7
module load samtools

#create directories
cd /pub/galentm/283/
mkdir bismarktest

#download maize reference genome and unzip
cd bismarktest/
mkdir ref
cd ref/
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
gunzip Zea_mays.AGPv4.dna.toplevel.fa.gz

#download B73 maize WGBS data
cd ../
mkdir reads
cd reads/
prefetch -v SRR408788
prefetch -v SRR408796
prefetch -v SRR408797
prefetch -v SRR408795
prefetch -v SRR408794
prefetch -v SRR408791
prefetch -v SRR408790
#downloaded to my home directory for some reason...

cd ~/ncbi/public/sra/
mv SRR4087* /pub/galentm/283/bismarktest/reads/

#Convert from .sra -> .fq
cd /pub/galentm/283/bismarktest/reads/
fastq-dump SRR4087*

#Trim reads (min length 30 bp according to Danelle's manuscript). 
#fastqc SRR408788.fastq <- outside script
#fastQC shows that the Illumina universal adapter is present in the reads
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408788.fastq SRR408788.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408796.fastq SRR408796.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408797.fastq SRR408797.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408794.fastq SRR408794.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408791.fastq SRR408791.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408790.fastq SRR408790.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20
java -jar /data/apps/trimmomatic/0.35/trimmomatic-0.35.jar SE -phred33 SRR408795.fastq SRR408795.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.35/adapters/TruSeq3-SE.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 AVGQUAL:20

#Genome indexing/bisulfite conversion. At this step, the reference genome is (computationally) bisulfite converted on both forward an reverse strands, creating two new reference genomes: one with C -> T converion and one with G -> A conversion
#This is done to all fasta files in directory
bismark_genome_preparation --verbose /pub/galentm/283/bismarktest/ref/

#Alignment with bismark/bowtie2
cd /pub/galentm/283/bismarktest/
bismark --bowtie2 --genome /pub/galentm/283/bismarktest/ref/ --se SRR408788.trimmed.fastq, SRR408796.trimmed.fastq, SRR408797.trimmed.fastq, SRR408795.trimmed.fastq, SRR408794.trimmed.fastq, SRR408791.trimmed.fastq, SRR408790.trimmed.fastq

#Extract methylation data for each cytosine position
bismark_methylation_extractor --report --multicore ${NSLOTS} --single-end --cytosine_report --CX --genome_folder /pub/galentm/283/bismarktest/ref/ *bt2.bam
