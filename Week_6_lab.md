
### Create projects
```
createProject ATACseq /pub/galentm/283/
createProject RNAseq /pub/galentm/283/
createProject DNAseq /pub/galentm/283/
```

### ATAC seq

````
cp * ../../DNAseq/data/raw/
cd ../../DNAseq/data/raw/


for file in  $(ls *.fq.gz)
do
  mv ${file} $(grep $(echo $file | cut -f5 -d"_") README.ATACseq.txt | awk '{print $2"_"$3"_"$4}')_$(echo $file | cut -f6 -d"_")
done
```

### RNA seq

```
cp Project_plex*/Sample_*/* ../../../RNAseq/data/raw/
cd ../
cp RNAseq384_SampleCoding.txt ../../RNAseq/data/raw/
cd ../../RNAseq/data/raw/

# No real need to rename anything here as far as I can tell
#for file in  $(ls *.fastq.gz)
#do
#  echo mv ${file} $(grep $(echo $file | cut -f1 -d"_") README.ATACseq.txt | awk '{print $2"_"$3"_"$4}')_#$(echo $file | cut -f6 -d"_")
#done
```

### DNA seq

```
cp * ../../DNAseq/data/raw/
cd ../../283/DNAseq/data/raw/


# Loops that rename arbitrary names with sample names

for file in $(ls ADL06*.fq.gz)
do
  mv ${file} ${file//ADL06/A4}
done

for file in $(ls ADL09*.fq.gz)
do
  mv ${file} ${file//ADL09/A5}
done

for file in $(ls ADL10*.fq.gz)
do
  mv ${file} ${file//ADL10/A6}
done

for file in $(ls ADL14*.fq.gz)
do
  mv ${file} ${file//ADL14/A7}
done
```

If you want to take a look in the project directories, they're all in /pub/galentm/283/ on the HPC
