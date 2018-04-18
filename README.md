PSI-Sigma Ψ<sub>Σ</sub>
=================
Previous single-exon PSI approaches were designed for simple splicing events with only one alternative exon, but they can be ambiguous in the case of mutually exclusive exons, multi-exon skipping, and more complex events. The PSI-sigma is using a new splicing index that is more flexible, can incoporate novel junctions, and can compute PSI values of individual exons in complex splicing events.

AUTHOR/SUPPORT
==============
Kuan-Ting (Woody) Lin, klin@cshl.edu

MANUAL
======
Step 0. Generate alignment files (.bam) by splice-aware alignment tools (e.g., STAR: https://github.com/alexdobin/STAR)
Step 1. Extract junction read information (If you don't have SJ.out.tab file from STAR)
```

mv <Junction Read Files> <Junction Folder>
```
Step 2. Build a database of splicing events for each chromosome (or download here:)
```
perl PSIsigma-db.pl <GTF File> <Junction Folder> <Chromosome>
cat chr*.db > <Database File>
```
Step 3. (Optional) Extract intronic read information
```
perl ir.pl <Database File> <BAM File> <Type>
mv <Intronic Read Files> <Junction Folder>
```
Step 4. Estimate PSI values for each splicing event
```
cd <Junction Folder>
perl PSIsigma-v.1.0.pl <Database File> <Output File>
```
Step 5. Filter and annotated splicing events
```
perl filter.pl GRCh38.87.mapping.txt <Output File>
```
PERFORMANCE
==============


SOFTWARE REQUIREMENTS
==============================
 * Perl (https://www.perl.org/get.html)

Perl EXTENTIONS
==============================
 * PDL::LiteF
 * PDL::Stats
 * PDL::GSL::CDF
 * Statistics::Multtest

INSTALLATION EXAMPLE
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
* Linux Bash Shell on Windows: https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/

LIMITATIONS
===========

CITATION
===========
Please cite: https://www.ncbi.nlm.nih.gov/pubmed/29449409
* Lin KT, Ma WK, Scharner J, Liu YR, Krainer AR. 2018. A human-specific switch of alternatively spliced AFMID isoforms contributes to TP53 mutations and tumor recurrence in hepatocellular carcinoma. Genome Res doi:10.1101/gr.227181.117.
