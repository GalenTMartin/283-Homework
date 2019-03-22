# Final Project: Aligning whole genome bisulfite sequencing reads and identifying methylated cytosines

### Rationale/background
This is a bit of a deviation from what I initially proposed to do (comparing expression data to CHH island-associated genes). Aside from the fact that this is more germane to what we've been learning in class, I decided to switch for two main reasons. Firstly, I realized that I'll be recieving whole genome bisulfite sequencing (WGBS) data from the inbred maize lines in a couple of months and I didn't have any idea how to process it. Secondly, as I've been working with Danelle's grass methylome data, I've found that my overall summary statistics aren't lining up with similar statistics from other papers which use the same data (for example, B73 maize should have ~5% CHH methylation, but it has only ~1.5% in Danelle's data. As a result, I'd like to reanalyze some of this data in different ways and see how different they are.

There are several programs which folks use to work with WGBS data. The predominant ones I found in papers were BS-seeker2, bismark, and MethylPy. I've chosen to start off with bismark since it's the program used by for the grass data that I've been looking at, but ultimately I'd like to reanalyze the data with a MethylPy pipeline as well since this was the program used by one of the main papers to which I've been comparing my data ([Niederhuth et al., 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1059-0)). 

### Steps
1. Download reads from NCBI database
2. Download reference genome
3. Bisulfite convert reference with bismark
4. Trim reads with trimmomatic
5. Align reads using bismark/bowtie2
6. Obtain methylated/nonmethylated reads data from bismark_methylation_extractor
7. Use binomial test to find real methylated cytosine positions


### Pipeline
The first part of the script downloads the files and creates directories for the project:
```

cd /pub/galentm/283/
mkdir bismarktest

cd bismarktest/
mkdir ref
cd ref/
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-36/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
gunzip Zea_mays.AGPv4.dna.toplevel.fa.gz
```

The second part retrieves the SRA files from NCBI. According to what I read online, it's possible to do this with fastq-dump <access #>, but this wasn't working for me, so I ended up using prefetch to download the regular .sra files. These files were downloaded to a new folder in my root directory, so I moved them to my project directory. These files were then converted from .sra to .fastq using fastq-dump. At this point I started just using accession SRR408788 in order to have enough time to finish.
```
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

cd ~/ncbi/public/sra/
mv SRR4087* /pub/galentm/283/bismarktest/reads/

fastq-dump SRR408788
```

Next, I trimmed the reads for quality/adapters. Initially I had trouble with Trimmomatic as a result of not reading the documentation carefully enough...
```
Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 1
        at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:54)
        at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:32)
        at org.usadellab.trimmomatic.Trimmomatic.createTrimmers(Trimmomatic.java:41)
        at org.usadellab.trimmomatic.TrimmomaticSE.run(TrimmomaticSE.java:298)
        at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:67)
```
This was fixed by adding the seed mismatches, palindrome clip threshold, and simple clip threshold to the ILLUMINACLIP option.

After trimming, I used a test dataset of the trimmed reads to test the bismark alignment function:
```
head -n 10000 SRR408788.trimmed2.fastq > testdata.trimmed.fastq
bismark --bowtie2 --genome /pub/galentm/283/bismarktest/ref/ --se testdata.trimmed.fastq
```
Bismark also prints a summary at this stage:
```
Final Cytosine Methylation Report
=================================
Total number of C's analysed:   802

Total methylated C's in CpG context:    206
Total methylated C's in CHG context:    138
Total methylated C's in CHH context:    27
Total methylated C's in Unknown context:        0

Total unmethylated C's in CpG context:  27
Total unmethylated C's in CHG context:  57
Total unmethylated C's in CHH context:  347
Total unmethylated C's in Unknown context:      0

C methylated in CpG context:    88.4%
C methylated in CHG context:    70.8%
C methylated in CHH context:    7.2%
```

Finally, I used bismark methylation extractor to get a file with the coverage and methylation state of each cytosine in the genome.
```
bismark_methylation_extractor --report --single-end --cytosine_report --CX --genome_folder /pub/galentm/283/bismarktest/ref/ testdata.trimmed_bismark_bt2.bam
```

The file generated looks like this:
```
6       76217380        -       1       0       CG      CGA
6       76217383        -       1       0       CG      CGT
6       76217386        -       1       0       CG      CGT
6       76217394        -       1       0       CG      CGA
6       76217395        -       1       0       CHG     CCG
6       76217399        -       1       0       CHG     CTG
6       76217414        -       1       0       CG      CGG
6       76217418        -       1       0       CHG     CAG
6       76217431        -       1       0       CG      CGC
6       76217437        -       1       0       CG      CGT
6       76217441        -       1       0       CG      CGC
6       89408276        +       1       0       CHG     CAG
6       89408285        +       1       0       CHG     CAG
6       89408306        +       1       0       CHG     CAG
6       94814140        -       1       0       CG      CGG
6       94814153        -       1       0       CHG     CAG
```
where the columns are chromosome, position, strand, count methylated, count unmethylated, context motif, and C context. 

The next step in the process will be to find the false methylation rate by mapping reads to the chloroplast genome (which is assumed to be unmethylated) and apply the false methylation rate through the binomial test used by [Lister et al. 2007](https://www.sciencedirect.com/science/article/pii/S0092867408004480) to include only cytosines which have a high probability of being methylated. 

