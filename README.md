PSI-Sigma
=================
Percent Spliced-In (PSI) values are commonly used to report alternative pre-mRNA splicing (AS) changes.
However, previous PSI-detection methods are limited to specific types of AS events. PSI-Sigma is using a new splicing index (PSI<sub>Σ</sub>) that is more flexible, can incoporate novel junctions, and can compute PSI values of individual exons in complex splicing events.
<br/>
 * PSI-Sigma is now released in obfuscated codes.

AUTHOR/SUPPORT
==============
Kuan-Ting (Woody) Lin, klin@cshl.edu

Alignment files
======
For short-read RNA-seq data, please generate .bam, .bai and .SJ.out files by using STAR (https://github.com/alexdobin/STAR).
```
###This is an example for short-read RNA-seq###
STAR --runThreadN 6 \
	--outSAMtype BAM SortedByCoordinate \
	--outFilterIntronMotifs RemoveNoncanonical \
	--genomeDir ~/index/starR100H38 \
	--twopassMode Basic \
	--readFilesIn R1.fastq R2.fastq \
	--outFileNamePrefix <NAME>.
samtools index <NAME>.Aligned.sortedByCoord.out.bam
```
For long-read RNA-seq data, please use GMAP (http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz).
```
###This is an example for long-read RNA-seq###
~/gmap-2017-11-15/bin/gmap -d GRCh38 -f samse --min-trimmed-coverage=0.5 --no-chimeras -B 5 -t 6 ~/MinION_long_read.fastq > <NAME>.sam
samtools view -bS <NAME>.sam > <NAME>.bam
samtools sort <NAME>.bam -o <NAME>.Aligned.sortedByCoord.out.bam
samtools index <NAME>.Aligned.sortedByCoord.out.bam
```
Quick Start
======
Create links to the .bam, bai, and .SJ.out files in the a folder (afolder). If you are using long-read RNA-seq data, .SJ.out files will be generated automatically since GMAP doesn't produce the file.
```
mkdir afolder
cd afolder
ln -s bamfolder/*.bam* .
ln -s bamfolder/*.SJ.* .
```
Download a .gtf file and sort the coordinates.
```
get ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens//Homo_sapiens.GRCh38.87.gtf.gz
gzip -d Homo_sapiens.GRCh38.87.gtf.gz
(grep "^#" Homo_sapiens.GRCh38.87.gtf; grep -v "^#" Homo_sapiens.GRCh38.87.gtf | sort -k1,1 -k4,4n) > Homo_sapiens.GRCh38.87.sorted.gtf
rm Homo_sapiens.GRCh38.87.gtf
```
Create two files: (1) groupa.txt and (2) groupb.txt. Please put suffixes of your .SJ.out.tab files in the groupa.txt or groupb.txt. For example, the suffix of a "Sequins_MixA.SJ.out.tab" file is "Sequins_MixA". Groupa.txt will be compared with groupb.txt. Below is an example of processing .SJ.out.tab files from TCGA (11A means normal and 01A means tumor):
```
#For TCGA files:
ls *-11A-*.SJ* | sed s/.SJ.out.tab//g > groupa.txt
ls *-01A-*.SJ* | sed s/.SJ.out.tab//g > groupb.txt

#Alternatively, you can just put the names of your .bam files:
echo Sequins_MixA.Aligned.sortedByCoord.out.bam > groupa.txt
echo Sequins_MixB.Aligned.sortedByCoord.out.bam > groupb.txt
```
Run dummyai.pl. After the .gtf file, please specify 1 for short-read RNA-seq and 2 for long-read RNA-seq. The last column is used to specify the minimum number of supporting reads for an AS event (10 is specified in the example below).
```
#For short-read RNA-seq (minimum 10 supporting reads for an AS event)
perl ~/PSIsigma/dummyai.pl Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 1 10
#For long-read RNA-seq (minimum 10 supporting reads for an AS event)
perl ~/PSIsigma/dummyai.pl Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 2 10
```
That's it.
The results will be in the PSIsigma_r3_ir10.sorted.txt.

 * Junction Read File: *.SJ.out.tab
 * Intronic Read File: *.IR.out.tab
 * Database File: *.db
 * BAM File: *.bam
 * GTF File: *.gtf (http://useast.ensembl.org/info/data/ftp/index.html/)


PERFORMANCE
==============


SOFTWARE REQUIREMENTS
==============================
 * Perl (https://www.perl.org/get.html)
 * Samtools (http://www.htslib.org)

Perl EXTENTIONS
==============================
 * PDL::LiteF
 * PDL::Stats
 * PDL::GSL::CDF
 * Statistics::Multtest

EXAMPLE of INSTALLING Perl EXTENTIONS
============================== 
```
# 0. Set up working directory for Perl library (Using Perl version 5.18 as an example)
export PERL5LIB=/usr/local/lib/perl/5.18

# 1. Install cpanm
cpan App::cpanminus
cpanm PDL::LiteF
cpanm PDL::Stats

# 2. Install GSL (Using GSL version 2.4 as an example)
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
tar zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure
make
make install
cd ..

# 3. Install PDL::GSL
cpanm PDL::GSL::CDF
cpanm Statistics::Multtest
```
PSI-Sigma on Windows OS
===========
PSI-Sigma has been tested in Linux and Mac OS environment. You can install Linux bash shell on Windows to run PSI-Sigma.
* Linux Bash Shell on Windows: https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/

LIMITATIONS
===========

CITATION
===========

PUBLICATIONS USING PSI-Sigma
===========
https://www.ncbi.nlm.nih.gov/pubmed/29449409
* Lin KT, Ma WK, Scharner J, Liu YR, Krainer AR. 2018. A human-specific switch of alternatively spliced AFMID isoforms contributes to TP53 mutations and tumor recurrence in hepatocellular carcinoma. Genome Res doi:10.1101/gr.227181.117.

Commercial Use
===========
* For licensing, please contact CSHL tech transfer office: narayan@cshl.edu

[![Analytics](https://ga-beacon.appspot.com/UA-123441271-1/PSI-Sigma/readme)](https://github.com/igrigorik/ga-beacon)

