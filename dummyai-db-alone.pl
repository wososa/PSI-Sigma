=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
	use strict;
	use Cwd qw(abs_path);
	
	my ($gtf,$name) = @ARGV;

	my $path = abs_path($0);
	$path=~s/\/dummyai\-db\-alone\.pl//;
	print "Path = $path\n";
	
	my $noveljunctioncriteria = 10;
	
	my %chr;
    open(FILE,"$gtf") || die "Aborting.. Can't open $gtf : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line=~/^\#/);
        my @array = split(/\t/,$line);
        next if($array[2] ne "transcript");
        my ($chr,$cat,$start,$end,$strand,$name) = ($array[0],$array[1],$array[3],$array[4],$array[6],$array[8]);
        $chr = "chr" . $chr if($chr!~/chr/);
        $chr{$chr}++;
    }
    close(OUT);
    close(FILE);
    
	my $starttime = time;
	my $dbname = $name . ".db";
	my $bedname = $name . ".bed";
	my $chrs;
	foreach my $chr(sort keys %chr){
		$chrs .= "\t" . $chr;
	}
	$chrs=~s/\t//;
	if(-e $dbname){
		if(-z $dbname){
			print "Regenerating $dbname...\n";
			rundb($noveljunctioncriteria,$gtf,$chrs);
		}
	}else{
		print "Generating $dbname...\n";
		rundb($noveljunctioncriteria,$gtf,$chrs);
	}
	my $stoptime = time;
    my $hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
    print "===Database spent $hours hours.===\n";

sub rundb{
	my $noveljunctioncriteria = shift;
	my $gtf = shift;
	my $chrs = shift;
	my @chromosomes = split(/\t/,$chrs);	
	foreach my $chromosome(@chromosomes){
		next if($chromosome=~/chrGL/);
		next if($chromosome=~/chrKI/);
		my $commend = "perl " . $path . "/PSIsigma-db-v.1.0.pl $gtf " . $chromosome . " " . $noveljunctioncriteria;
		#print "Doing... $commend\n";
		print "Doing... $chromosome\n";
		system("$commend");
	}
	
	system("cat chr*.db > $dbname");
	system("cat chr*.bed > $bedname");
	system("rm chr*.db");
	system("rm chr*.bed");
}
	