=begin
PSI-Sigma
A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2022

PSI-Sigma is a free open source software distributed under GPLv3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
=end
=cut
#!/usr/bin/perl -w
	use strict;
	use Cwd qw(abs_path);
	
	my ($gtf,$name,$longread,$supporting_read_criteria) = @ARGV;
	
	my $path = abs_path($0);
	$path=~s/\/dummyai\.pl//;
	print "Path = $path\n";
	    
    my ($starttime,$stoptime,$hours) = (time,0,0);
    
	if(scalar @ARGV != 4){
		print "Please specify .gtf file, output name, long(2) or short(1) read, and number of supporting reads.\n";
		exit;
	}

	my $noveljunctioncriteria = 10;
	#my $supporting_read_criteria = 10;
	my $intron_criteria = 3;
	my $totaltime = 0;
	
    my %group;
 	open(FILE,"groupa.txt") || die "Aborting.. Can't open groupa.txt : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        $group{$line}++;
    }
 	open(FILE,"groupb.txt") || die "Aborting.. Can't open groupb.txt : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        $group{$line}++;
    }

	print "Generating mapping file...\n";
    my $mappingcount = 0;
	my %chr;
    open(FILE,"$gtf") || die "Aborting.. Can't open $gtf : $!\n";
    open(OUT,">" . $gtf . ".mapping.txt") || die "Aborting.. Can't open $gtf : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line=~/^\#/);
        my @array = split(/\t/,$line);
        next if($array[2] ne "transcript");
        my ($chr,$cat,$start,$end,$strand,$name) = ($array[0],$array[1],$array[3],$array[4],$array[6],$array[8]);
        $chr = "chr" . $chr if($chr!~/chr/);
        $chr{$chr}++;
		my $ENST = $name;
 		$ENST=~s/(.*)transcript\_id \"//;
    	$ENST=~s/\"\; (.*)//;
		if($ENST=~/\_/){
			$ENST=~s/\_/\./g;
		}
        $name=~s/(.*)gene\_name \"//;
    	$name=~s/\"\; (.*)//;
    	$name=~s/gene\_id \"//;
        print OUT $ENST . "\t" . $name . "\t" . $strand . "\n";
        $mappingcount++ if($ENST ne "" && $name ne "" && $strand ne "");
    }
    close(OUT);
    close(FILE);
    if($mappingcount == 0){
    	print "$gtf is not in an acceptable gtf format. Exiting...\n";
    	exit;
    }
    
    my $nfiles = scalar keys %group;
    if($nfiles < 2){
    	print "Not enough files. Exiting...\n";
		exit;
	}

	print "Checking splice-junction files...\n";
	$starttime = time;
	my $sjcount = 0;
	foreach my $bam(keys %group){
		next if($bam eq "");
		my $accession = $bam;
		$accession=~s/Aligned\.sortedByCoord\.out\.bam//;
		$accession=~s/sorted\.out\.bam//;
		$accession=~s/\.bam//;
		$accession=~s/\.$//;
		my $sjfn = $accession . ".SJ.out.tab";
		my $commend = "samtools view $bam | awk -f " . $path . "/sjFromSAMcollapseUandM_inclOverlaps.awk > " . $sjfn;
		if(-e $sjfn){
			if(-z $sjfn){
				print "Regenerating .SJ.out for $accession\n";
				if($bam=~ /Aligned\.sortedByCoord\.out\.bam/){
					print "$bam looks sorted by coordiates.\n";
					print "sjFromSAMcollapseUandM_inclOverlaps.awk may not generate accurate .SJ.out files. Please use .SJ.out from STAR aligner or sort the file by:\n";
					print "samtools sort -n -o $accession.Aligned.sortedByName.out.bam $bam\n";
					print "=CAUTION=[.IR.out will need a bam file sorted by coordiates]\n";
					#next;
				}
				system("$commend");
			}
		}else{
			print "Generating .SJ.out for $accession\n";
			if($bam=~ /Aligned\.sortedByCoord\.out\.bam/){
				print "$bam looks sorted by coordiates.\n";
				print "sjFromSAMcollapseUandM_inclOverlaps.awk may not generate accurate .SJ.out files. Please use .SJ.out from STAR aligner or sort the file by:\n";
				print "samtools sort -n -o $accession.Aligned.sortedByName.out.bam $bam\n";
				print "=CAUTION=[.IR.out will need a bam file sorted by coordiates]\n";
				#next;
			}
			system("$commend");
		}
		$sjcount++;
	}
	$stoptime = time;
    $hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
    $totaltime += $hours;
    print "===Splice-junction files spent $hours hours.===\n";
	if($sjcount < 2){
		print "Only $sjcount samples. It is not enough. Exiting...\n";
		exit;
	}
	
	$starttime = time;
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
	$stoptime = time;
    $hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
    $totaltime += $hours;
    print "===Database spent $hours hours.===\n";
    
    print "Getting intron reads....\n";
    $starttime = time;
    my $ircount = 0;
    foreach my $bam(keys %group){
		my $accession = $bam;
		$accession=~s/Aligned\.sortedByCoord\.out\.bam//;
		$accession=~s/sorted\.out\.bam//;
		$accession=~s/\.bam//;
		$accession=~s/\.$//;
		my $pname = "";
		my $checkpname = 0;
		if($bam!~/\.bam/){
			$pname = $bam . "\.Aligned\.sortedByCoord\.out\.bam";
			if(-e $pname){
				$bam = $pname;
				$checkpname = 1;
			}
			if($checkpname == 0){
				$pname = $bam . "\.sorted\.out\.bam";
				if(-e $pname){
					$bam = $pname;
					$checkpname = 1;
				}
			}
			if($checkpname == 0){
				$pname = $bam . "\.bam";
				if(-e $pname){
					$bam = $pname;
					$checkpname = 1;
				}
			}
		}
		print "Checking $bam...\n";
		my $irfn = $accession . ".IR.out.tab";
		print "Checking $irfn...\n";
		my $commend = "perl " . $path . "/PSIsigma-ir-v.1.0.pl " . $name . ".db " . $bam . " " . $longread;
		#print "commend = $commend\n";
		if(-e $irfn){
			if(-z $irfn){
				print "Regenerating .IR.out for $accession\n";
				system("$commend");
			}else{
				print "$irfn existed. Pass...\n";
			}
		}else{
			print "Generating .IR.out for $accession\n";
			system("$commend");
		}
		$ircount++;
    }
	$stoptime = time;
	$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
	$totaltime += $hours;
	print "===Intron-read file spent $hours hours.===\n";
	if($ircount < 2){
		print "Only $ircount samples. It is not enough. Exiting...\n";
		exit;
	}
	
	print "Ready to do PSI analysis...\n";
	$starttime = time;
	my $commend = "perl " . $path . "/PSIsigma-PSI-v.1.0.pl " . $name . ".db " . $name . " " . $supporting_read_criteria . " " . $intron_criteria . " " . $longread;
	system($commend);
	$stoptime = time;
	$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
	$totaltime += $hours;
	print "===PSI analysis spent $hours hours.===\n";
	
	print "Filtering ΔPSI results...\n";
	$starttime = time;
	my $nofilter = 0;
	$nofilter = 1 if($nfiles < 4);
	$commend = "perl " . $path . "/PSIsigma-filter-v.1.0.pl " . $name . ".db " . $gtf . ".mapping.txt " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt " . $nofilter;
	system($commend);
	$stoptime = time;
	$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
	$totaltime += $hours;
	print "===Filtering spent $hours hours.===\n";
	
	print "\n***Total: $totaltime hours (or " . ($totaltime*60) . "mins).\n";
	
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
	