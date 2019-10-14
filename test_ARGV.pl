=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
	use strict;
	use PDL::LiteF;
	use PDL::Stats;
	use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
	
    my ($db,$outputassccession,$criteria,$introallcriteria,$longread) = @ARGV;
    
    if(scalar @ARGV != 5){
    	print "Please specify 5 parameters: (1) database, (2) output name , (3) minimum supporting junction reads, (4) minimum intron supporting reads, and (5) if the input data is short or long read.\n"; 
    	exit;
    }
	print "$db,$outputassccession,$criteria,$introallcriteria,$longread\n";