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
#one file name per line in groupa.txt and groupb.txt
echo A1.sortedByCoord.out.bam > groupa.txt
echo B1.sortedByCoord.out.bam > groupb.txt
```
Ready to go
======
You need to generate PSIsigma.db first:
```
perl ~/PSIsigma/dummyai-db-alone.pl Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma
```
Next, you can generate IR.out.tab files for each .bam files, respectively:
```
perl ~/PSIsigma/PSIsigma-ir-v.1.0.pl PSIsigma.db A1.sortedByCoord.out.bam 1
perl ~/PSIsigma/PSIsigma-ir-v.1.0.pl PSIsigma.db B1.sortedByCoord.out.bam 1
<repeat for all .bam files in groupa.txt and groupb.txt>
```
After obtaining the .db and .IR.out.tab files, you can run dummyai.pl.
```
#For short-read RNA-seq (minimum 10 supporting reads for an AS event)
perl ~/PSIsigma/dummyai.pl Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 1 10
#For long-read RNA-seq (minimum 10 supporting reads for an AS event)
perl ~/PSIsigma/dummyai.pl Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 2 10
```
That's it. The results will be in the PSIsigma_r10_ir3.sorted.txt.

