## Alignments and indexing 

### Indexing: 
```
ref="/pub/galentm/283/dmel-all-chromosome-r6.13.fasta"
bwa index $ref
samtools faidx $ref
java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=$ref O=/pub/galentm/283/dmel-all-chromosome-r6.13.fasta.dict
bowtie2-build $ref
```

### DNAseq alignment
```
#!/usr/bin/env bash
#
#$ -N w7DNA_1
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
module load samtools/1.3
module load bcftools/1.3
module load enthought_python/7.3.2
module load gatk/2.4-7
module load picard-tools/1.87
module load java/1.7
module load bowtie2

ref="/pub/galentm/283/dmel-all-chromosome-r6.13.fasta"
file="/pub/galentm/283/DNAseq/data/processed/fixed2*fq.gz"
folder="/pub/galentm/283/DNAseq/data/processed/"

cd /pub/galentm/283/DNAseq/data/processed

#bwa index ${ref}
#samtools faidx  ${ref}
#java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=${ref} O=/pub/galentm/283/dmel-all-chromosome-r6.13.fasta.dict
#bowtie2-build ${ref}


#ls *1.fq.gz | sed 's/_1.fq.gz//' >DNAseq.prefixes.txt  # outside the array
prefix=`head -n ${SGE_TASK_ID} DNAseq.prefixes.txt | tail -n 1`
bwa mem -t 8 -M ${ref} ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -bS - > ${prefix}.bam

echo $SGE_TASK_ID > SGEtest.txt
#ls *1.fq.gz | sed 's/_1.fq.gz//' > DNAseq.prefixes.txt
#prefix=`head -n $SGE_TASK_ID DNAseq.prefixes.txt | tail -n 1`
#bwa mem -t 8 -M $ref $prefix_1.fq.gz $prefix_2.fq.gz | samtools view -bS - > folder/$SGE_TASK_ID.bam
#samtools sort folder/$SGE_TASK_ID.bam -o folder/$SGE_TASK_ID.sort.bam

java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=${prefix}.sort.bam O=${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT

samtools index ${prefix}.RG.bam

echo ${ref} > reftest.txt
echo ${prefix} >> prefixtest.txt
```

### RNAseq alignment
```
#!/usr/bin/env bash
#
#$ -N w7RNA
#$ -t 1
#$ -R y
#$ -q epyc,free64,abio,bio,bsg,bsg2
#$ -m beas
#$ -M galentm@uci.edu
### -o
### -e
#$ -pe openmp 4
#$ -ckpt restart

module load bwa/0.7.8
module load bcftools/1.3
module load picard-tools/1.87
module load gatk/2.4-7
module load java/1.7
module load samtools/1.3
module load tophat/2.1.0
module load bowtie2/2.2.7

#ls *.fastq.gz | sed 's/_R1_001.fastq.gz//' > /pub/galentm/283/ATACseq/data/processed/RNAseq.prefixes.txt

prefix=`head -n ${SGE_TASK_ID} RNAseq.prefixes.txt | tail -n 1`

tophat -p 8 -G /pub/galentm/283/dmel-all-r6.13.gtf -o /pub/galentm/283/RNAseq/data/processed /pub/galentm/283/dmel-all-chromosome-r6.13.fasta.btout ${prefix}_R1_0001.fastq.gz ${prefix}_R2_001.fastq.gz

samtools sort /pub/galentm/283/RNAseq/data/processed/accepted_hits.bam -o /pub/galentm/283/RNAseq/data/processed/accepted_hits.sort.bam
samtools index
```

### ATACseq alignment
```
#!/usr/bin/env bash
#
#$ -N align030219
#$ -t 1-24
#$ -R y
#$ -q epyc,free64,abio,bio,bsg,bsg2
#$ -m beas
#$ -M galentm@uci.edu
### -o
### -e
#$ -pe openmp 32-128
#$ -ckpt restart

module load bwa/0.7.8
module load bcftools/1.3
module load picard-tools/1.87
module load gatk/2.4-7
module load java/1.7
module load samtools/1.3

cd /pub/galentm/283/ATACseq/data/processed

ref="/pub/galentm/283/dmel-all-chromosome-r6.13.fasta"
file="/pub/galentm/283/ATACseq/data/processed/fixed2*fq.gz"
folder="/pub/galentm/283/ATACseq/data/processed/"

#bwa index ${ref}
#samtools faidx  ${ref}
#java -d64 -Xmx128g -jar /data/apps/picard-tools/1.87/CreateSequenceDictionary.jar R=${ref} O=/pub/galentm/283/dmel-all-chromosome-r6.13.fasta.dict
#bowtie2-build ${ref}


#ls *1.fq.gz | sed 's/_1.fq.gz//' >DNAseq.prefixes.txt  # outside the array
prefix=`head -n ${SGE_TASK_ID} ATACseq.prefixes.txt | tail -n 1`
bwa mem -t 8 -M ${ref} ${prefix}_R1.fq.gz ${prefix}_R2.fq.gz | samtools view -bS - > ${prefix}.bam

echo $SGE_TASK_ID >> SGEtest.txt

samtools sort ${prefix}.bam -o ${prefix}.sort.bam

java -Xmx20g -jar /data/apps/picard-tools/1.87/AddOrReplaceReadGroups.jar I=${prefix}.sort.bam O=${prefix}.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${prefix} RGSM=${prefix} VALIDATION_STRINGENCY=LENIENT

samtools index ${prefix}.RG.bam

#ls *1.fq.gz | sed 's/_1.fq.gz//' > /pub/galentm/283/ATACseq/data/processed/ATACseq.prefixes.txt
```

