
#!/usr/bin/env bash
#
#$ -N labweek7_1
#$ -t 1-12
#$ -R y
#$ -q epyc,free64,abio,bio
#$ -m beas
#$ -M galentm@uci.edu
### -o
### -e
#$ -pe openmp 2
#$ -ckpt restart

module load bwa/0.7.8
module load bcftools/1.3
module load picard-tools/1.87
module load gatk/2.4-7
module load java/1.7
module load samtools/1.3

cd /pub/galentm/283/DNAseq/data/processed

ref="/pub/galentm/283/dmel-all-chromosome-r6.13.fasta"
file="/pub/galentm/283/DNAseq/data/processed/fixed2*fq.gz"
folder="/pub/galentm/283/DNAseq/data/processed/"

bwa index ${ref}
samtools faidx  ${ref}
java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=${ref} O=/pub/galentm/283/dmel-all-chromosome-r6.13.fasta.dict
bowtie2-build ${ref}


#ls *1.fq.gz | sed 's/_1.fq.gz//' >DNAseq.prefixes.txt  # outside the array
prefix=`head -n ${SGE_TASK_ID} DNAseq.prefixes.txt | tail -n 1`
bwa mem -t 8 -M ${ref} ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -bS - > ${prefix}.bam

echo $SGE_TASK_ID > SGEtest.txt
#ls *1.fq.gz | sed 's/_1.fq.gz//' > DNAseq.prefixes.txt
#prefix=`head -n $SGE_TASK_ID DNAseq.prefixes.txt | tail -n 1`
#bwa mem -t 8 -M $ref $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -bS - > folder/$SGE_TASK_ID.bam
#samtools sort folder/$SGE_TASK_ID.bam -o folder/$SGE_TASK_ID.sort.bam
#java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=folder/$SGE_TASK_ID.sort.bam O=folder/$SGE_TASK_ID.RG.bam SORT_ORDER=coordinate #RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID={your sample name} RGSM={your sample name} VALIDATION_STRINGENCY=LENIENT
#samtools index folder/$prefix.RG.bam

samtools sort ${prefix}.bam -o ${prefix}.sort.bam

java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=${prefix}.sort.bam O=${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT

samtools index ${prefix}.RG.bam

echo ${ref} > reftest.txt
echo ${prefix} >> prefixtest.txt
