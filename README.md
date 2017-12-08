The New PSI index
=================
Previous single-exon PSI approaches were designed for simple splicing events with only one alternative exon, but they can be ambiguous in the case of mutually exclusive exons, multi-exon skipping, and more complex events. The new PSI index is flexible, can incoporate novel junctions, and can compute PSI values of individual exons in complex splicing events.

AUTHOR/SUPPORT
==============
Kuan-Ting (Woody) Lin, klin@cshl.edu

MANUAL
======

PERFORMANCE
==============


SOFTWARE REQUIREMENTS
==============================
  * Alignment files (.bam) from splice-aware alignment tools (e.g., STAR: https://github.com/alexdobin/STAR)
  * Perl (https://www.perl.org/get.html)

Perl EXTENTIONS
==============================
  * PDL::LiteF
  * PDL::Stats

Linux, Mac and "Linux Bash Shell on Windows"
-----

```
# Install Perl extentions
cpan App::cpanminus
cpanm PDL::LiteF
cpanm PDL::Stats

#Alternatively,
cpan
o conf urllist ftp://cpan.hexten.net/ ftp://mirrors.rit.edu/CPAN/ http://mirror.nyi.net/CPAN/
o conf commit
install PDL::LiteF
install PDL::Stats

# If the urllist is slow for you, please visit http://www.cpan.org/SITES.html and pick the ones closer to your region.
# Linux Bash Shell on Windows: https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/

```


LIMITATIONS
===========



