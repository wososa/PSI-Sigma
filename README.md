PSI<sub>Σ</sub>
=================
Previous single-exon PSI approaches were designed for simple splicing events with only one alternative exon, but they can be ambiguous in the case of mutually exclusive exons, multi-exon skipping, and more complex events. PSI-Sigma is using a new splicing index (Ψ<sub>Σ</sub>) that is more flexible, can incoporate novel junctions, and can compute PSI values of individual exons in complex splicing events.

AUTHOR/SUPPORT
==============
Kuan-Ting (Woody) Lin, klin@cshl.edu

MANUAL
======
Generate .bam files by splice-aware alignment tools (e.g., STAR: https://github.com/alexdobin/STAR) and create a folder (afolder) for the links to the .bam files.
```
mkdir afolder
cd afolder
ln -s bamfolder/*.bam .
```
Download a .gtf file and sort the coordinates.
```
get ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens//Homo_sapiens.GRCh38.87.gtf.gz
gzip -d Homo_sapiens.GRCh38.87.gtf.gz
(grep "^#" Homo_sapiens.GRCh38.87.gtf; grep -v "^#" Homo_sapiens.GRCh38.87.gtf | sort -k1,1 -k4,4n) > Homo_sapiens.GRCh38.87.sorted.gtf
rm Homo_sapiens.GRCh38.87.gtf
```
Create two files: (1) groupa.txt and (2) groupb.txt. The .bam files in groupa.txt will be compared with groupb.txt.
The examples here are assuming for .bam files downloaded from TCGA:
```
ls *-11A-*.bam > groupa.txt
ls *-01A-*.bam > groupb.txt
```
Run dumpai.pl and specify the folder (e.g., ~/PSIsigma) where you put the PSIsigma scripts.
Please specify 1 for short-read RNA-seq and 2 for long-read RNA-seq:
```
#For short-read RNA-seq
perl ~/PSIsigma/dumpai.pl ~/PSIsigma Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 1
#For long-read RNA-seq
perl ~/PSIsigma/dumpai.pl ~/PSIsigma Homo_sapiens.GRCh38.87.sorted.gtf PSIsigma 2
```
That's it.
The results will be in the PSIsigma_r10_ir3.filtered.txt.

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
 * threads

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
cpanm threads
```
* Linux Bash Shell on Windows: https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/

LIMITATIONS
===========

CITATION
===========
Please cite: https://www.ncbi.nlm.nih.gov/pubmed/29449409
* Lin KT, Ma WK, Scharner J, Liu YR, Krainer AR. 2018. A human-specific switch of alternatively spliced AFMID isoforms contributes to TP53 mutations and tumor recurrence in hepatocellular carcinoma. Genome Res doi:10.1101/gr.227181.117.

Commercial Uses
===========
Multi-thread pipeline: smartai.pl
Excel export
Confidence score
