Generating IR.out.tab file in parallel
=================
PSI-Sigma is currently a single-thread software. The majority of the computing time was used to create IR.out.tab files (intronic read counts). To generate IR.out.tab files in parallel, here are the instructions:

Before generating IR.out.tab files
======
Please complete the Quick Start steps (paste below):
```
mkdir afolder
cd afolder
ln -s bamfolder/*.bam* .
ln -s bamfolder/*.SJ.* .
get ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens//Homo_sapiens.GRCh38.87.gtf.gz
gzip -d Homo_sapiens.GRCh38.87.gtf.gz
(grep "^#" Homo_sapiens.GRCh38.87.gtf; grep -v "^#" Homo_sapiens.GRCh38.87.gtf | sort -k1,1 -k4,4n) > Homo_sapiens.GRCh38.87.sorted.gtf
rm Homo_sapiens.GRCh38.87.gtf
#For TCGA files:
ls *-11A-*.SJ* | sed s/.SJ.out.tab//g > groupa.txt
ls *-01A-*.SJ* | sed s/.SJ.out.tab//g > groupb.txt
#Alternatively, you can just put the names of your .bam files:
echo Sequins_MixA.Aligned.sortedByCoord.out.bam > groupa.txt
echo Sequins_MixB.Aligned.sortedByCoord.out.bam > groupb.txt
```
Ready to go
======
You need to generate PSIsigma.db first:
```

```
Next, you can generate IR.out.tab files for each .bam files, respectively:
```

```
That's it.

