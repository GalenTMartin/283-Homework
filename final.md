## Final Project: Aligning whole genome bisulfite sequencing reads and identifying methylated cytosines

### Rationale/background
This is a bit of a deviation from what I initially proposed to do (comparing expression data to CHH island-associated genes). Aside from the fact that this is more germane to what we've been learning in class, I decided to switch for two main reasons. Firstly, I realized that I'll be recieving whole genome bisulfite sequencing (WGBS) data from the inbred maize lines in a couple of months and I didn't have any idea how to process it. Secondly, as I've been working with Danelle's grass methylome data, I've found that my overall summary statistics aren't lining up with similar statistics from other papers which use the same data (for example, B73 maize should have ~5% CHH methylation, but it has only ~1.5% in Danelle's data. As a result, I'd like to reanalyze some of this data in different ways and see how different they are.

There are several programs which folks use to work with WGBS data. The predominant ones I found in papers were BS-seeker2, bismark, and MethylPy. I've chosen to start off with bismark since it's the program used by for the grass data that I've been looking at, but ultimately I'd like to reanalyze the data with a MethylPy pipeline as well since this was the program used by one of the main papers to which I've been comparing my data (Niederhuth et al., 2016). 

### Steps
1. Download reads from NCBI database
2. Download reference genome
3. Bisulfite convert reference with bismark
4. Trim reads with trimmomatic
5. Align reads using bismark/bowtie2
6. Obtain methylated/nonmethylated reads data from bismark _methylation_extractor
7. Use binomial test to find real methylated cytosine positions


### Analysis
Unfortunately, my bisulfite conversion of the reference still hasn't finished, though I've gotten it to work with the [test data set](https://github.com/FelixKrueger/Bismark/blob/master/test_data.fastq) given by the bismark authors. Currently I'm also stuck on getting timmomatic to work since it's giving me a java error:
```
Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 1
        at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:54)
        at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:32)
        at org.usadellab.trimmomatic.Trimmomatic.createTrimmers(Trimmomatic.java:41)
        at org.usadellab.trimmomatic.TrimmomaticSE.run(TrimmomaticSE.java:298)
        at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:67)
```
I realize that this has something to do with the fasta I'm giving it with adapter sequences, but I haven't had the time to work through it.
