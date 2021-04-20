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
	
	my $status = param(\@ARGV);
	
	my @param = split(/\t/,$status);
	
	if($param[0] eq "Help" || scalar @param < 4){
		print "\n";
		print "Usage: perl dummyai.pl [parameters]\n";
		print "Parameters:\n";
		print "  --gtf [text]		the gene annotation file for building PSI-Sigma database.\n";
		print "  --name [text]		the prefix of PSI-Sigma database and output files.\n";
		print "  --type [number]	1: short-read RNA-seq data\n";
		print "			2: long-read RNA-seq data\n";
		print "  --nread [number]	the minimum number of supporting reads for a splicing event.\n";
		print "  --skipratio [number]	the ratio (0~1) of skipping reads in a exon-skipping event [default: 0.05]. A splicing event has to have at least one input sample passing this criteria to be included in the analysis. \n";
		print "  --fmode [number]	0: delta-PSI > 10% and p-value < 0.01 (default/recommended)\n";
		print "			1: delta-PSI > 10%\n";
		print "			2: p-value < 0.05\n";
		print "			3: report all events\n";
		print "  --irmode [number]	0: only introns in the .gtf. (default)\n";
		print "			1: aggressively search for all introns\n";
		print "			2: Skip intron-retention events.\n";
		print "  --adjp [number]	0: Skip p-value adjustment. (default)\n";
		print "			1: Benjamini-Hochberg (Statistics::Multtest Perl module)\n";
		print "			2: Benjamini-Hochberg (qvalue or p.adjust() R package)\n";
		print "  --trimp [number]	Set FDR to 'na' for events whose maximum PSI values are below [number]% or whose minimum PSI values are above [100-number]% in all groupa.txt and groupb.txt during p-value adjustment. [default:5]\n";
		print "  --denominator [number]	0: Don't report denominators. (default)\n";
		print "			1: Report the table of denominators.\n";
		print "  --irrange [number]	0: IR event is using reference points of the target exon region.\n";
		print "			n: Use +-n bases around target intron region and +-n bases around event region. (default)\n";
		print "\n";
		exit;
	}
	
	my ($gtf,$name,$type,$supporting_read_criteria,$fmode,$skipratio,$irmode,$adjp,$trimp,$denominator,$irrange) = split(/\t/,$status);
	$fmode = 0 if($fmode ne "0" && $fmode ne "1" && $fmode ne "2" && $fmode ne "3");
	$irmode = 0 if($irmode ne "0" && $irmode ne "1" && $irmode ne "2");
	$irrange = 5 if(!$irrange || $irrange!~/\d/ || $irrange eq "-");
	$adjp = 0 if(!$adjp);
	$adjp = 0 if($adjp ne "0" && $adjp ne "1" && $adjp ne "2");
	$denominator = 0 if($denominator ne "0" && $denominator ne "1");
	$skipratio = 0.05 if($skipratio eq "-" || $skipratio > 1 || $skipratio < 0);
	$trimp = 5 if($trimp eq "-");
	
	print "gtf = $gtf\n";
	print "name = $name\n";
	print "type = $type\n";
	print "nread = $supporting_read_criteria\n";
	print "skipratio = $skipratio\n";
	print "fmode = $fmode\n";
	print "irmode = $irmode\n";
	print "adjp = $adjp\n";
	print "trimp = $trimp\n";
	print "denominator = $denominator\n";
	print "irrange = $irrange\n";
	
	if($adjp == 2){
		print "### NOTE: Staitistics::R and R are required.\n";
		use Statistics::R;
	}
	#my ($gtf,$name,$longread,$supporting_read_criteria) = @ARGV;
	
	my $path = abs_path($0);
	$path=~s/\/dummyai\.pl//;
	print "Path = $path\n";
	
    my ($starttime,$stoptime,$hours) = (time,0,0);

	my $noveljunctioncriteria = 10;
	#my $supporting_read_criteria = 10;
	my $intron_criteria = 3;
	my $totaltime = 0;
	
    my %group;
    my $patha = 0;
 	open(FILE,"groupa.txt") || die "Aborting.. Can't open groupa.txt : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        my $bam = $line;
        my $bamfn = $bam;
        my $bai = "$bam\.bai";
        my $tmp = "-";
        if($bam!~/\.bam/){
        	$tmp = $bam . ".Aligned.sortedByCoord.out.bam";
        	$bam = $tmp if(-e $tmp);
        	$tmp = $bam . ".sorted.bam";
        	$bam = $tmp if(-e $tmp);
        	$tmp = $bam . ".bam";
        	$bam = $tmp if(-e $tmp);
        }
        if($bam=~/\//){
        	$patha = 1;
        	$bamfn=~s/(.*)\///;
        	system("ln -s $bam $bamfn");
        	system("ln -s $bam\.bai $bamfn\.bai");
        	$bam = $bamfn;
        }
        if($bam=~/\.bam$/){
        	$bai = "$bam\.bai";
        	if(-e $bai){
        		print "$bai is ready.\n";
        		#system("ln -s $bam\.bai $bamfn\.bai");
        	}else{
        		print "$bai doesn't exist. Creating a new index...\n";
        		system("samtools index $bamfn");
        	}
        }
        $group{$bam}++;
    }
    
 	open(FILE,"groupb.txt") || die "Aborting.. Can't open groupb.txt : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        my $bam = $line;
        my $bamfn = $bam;
        my $bai = "$bam\.bai";
        my $tmp = "-";
        if($bam!~/\.bam/){
        	$tmp = $bam . ".Aligned.sortedByCoord.out.bam";
        	$bam = $tmp if(-e $tmp);
        	$tmp = $bam . ".sorted.bam";
        	$bam = $tmp if(-e $tmp);
        	$tmp = $bam . ".bam";
        	$bam = $tmp if(-e $tmp);
        }
        if($bam=~/\//){
        	$patha = 1;
        	$bamfn=~s/(.*)\///;
        	system("ln -s $bam $bamfn");
        	system("ln -s $bam\.bai $bamfn\.bai");
        	$bam = $bamfn;
        }
        if($bam=~/\.bam$/){
        	$bai = "$bam\.bai";
        	if(-e $bai){
        		print "$bai is ready.\n";
        		#system("ln -s $bam\.bai $bamfn\.bai");
        	}else{
        		print "$bai doesn't exist. Creating a new index...\n";
        		system("samtools index $bamfn");
        	}
        }
        $group{$bam}++;
    }
    
    if($patha == 1){
    	print "Formatting groupa.txt and groupb.txt...\n";
    	system("cp groupa.txt groupa.txt.orig");
    	system("cp groupb.txt groupb.txt.orig");
    	open(FILE,"groupa.txt.orig") || die "Aborting.. Can't open groupa.txt.orig : $!\n";
    	open(OUT,">groupa.txt") || die "Aborting.. Can't open groupa.txt : $!\n";
    	while(my $line=<FILE>){
        	chomp $line;
        	next if($line eq "");
        	$line=~s/(.*)\///;
        	$line=~s/\.Aligned\.sortedByCoord\.out\.bam//;
        	$line=~s/\.sorted\.out\.bam//;
        	$line=~s/\.sorted\.bam//;
        	$line=~s/\.bam//;
        	print OUT "$line\n";
        }
        close(OUT);
        close(FILE);
    	open(FILE,"groupb.txt.orig") || die "Aborting.. Can't open groupb.txt.orig : $!\n";
    	open(OUT,">groupb.txt") || die "Aborting.. Can't open groupb.txt : $!\n";
    	while(my $line=<FILE>){
        	chomp $line;
        	next if($line eq "");
        	$line=~s/(.*)\///;
        	$line=~s/\.Aligned\.sortedByCoord\.out\.bam//;
        	$line=~s/\.sorted\.out\.bam//;
        	$line=~s/\.sorted\.bam//;
        	$line=~s/\.bam//;
        	print OUT "$line\n";
        }
        close(OUT);
        close(FILE);
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
        #$chr = "chr" . $chr if($chr!~/chr/);
        $chr{$chr}++;
		my $ENST = $name;
 		$ENST=~s/(.*)transcript\_id \"//;
    	$ENST=~s/\"\; (.*)//;
		if($ENST=~/\_/){
			$ENST=~s/\_/\./g;
		}
		$name=~s/(.*)gene\_name \"//;
        $name=~s/\"\; (.*)//;
        $name=~s/\"\;//;
        $name=~s/gene\_id \"//;
        if($name=~/\_/){
			$name=~s/\_/\./g;
		}
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
        $accession=~s/\.Aligned\.sortedByCoord\.out\.bam//;
        $accession=~s/\.sorted\.out\.bam//;
        $accession=~s/\.sorted\.bam//;
		$accession=~s/\.bam//;
		$accession=~s/\.$//;
		my $sjfn = $accession . ".SJ.out.tab";
		
		my $commend = "samtools view $bam | awk -f " . $path . "/sjFromSAMcollapseUandM_inclOverlaps.awk > " . $sjfn;
		if(-e $sjfn){
			if(-z $sjfn){
				if($type == 1){
					generateSJ($bam,$accession);
				}
				if($type == 2){
					print "Generating... $sjfn\n";
					system($commend);
				}
			}
		}else{
			if($type == 1){
				generateSJ($bam,$accession);
			}
			if($type == 2){
				print "Generating... $sjfn\n";
				system($commend);
			}
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
			rundb($noveljunctioncriteria,$gtf,$chrs,$type);
		}
	}else{
		print "Generating $dbname...\n";
		rundb($noveljunctioncriteria,$gtf,$chrs,$type);
	}
	$stoptime = time;
    $hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
    $totaltime += $hours;
    print "===Database spent $hours hours.===\n";
    
    my $ircheck = 1;
    $ircheck = 0 if($irmode == 2);
    if($ircheck == 1){
		print "Getting intron reads....\n";
		$starttime = time;
		my $ircount = 0;
		foreach my $bam(keys %group){
			my $accession = $bam;
			$accession=~s/\.Aligned\.sortedByCoord\.out\.bam//;
			$accession=~s/\.sorted\.out\.bam//;
			$accession=~s/\.sorted\.bam//;
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
			my $commend = "perl " . $path . "/PSIsigma-ir-v.1.2.pl " . $name . ".db " . $bam . " " . $type;
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
	}
	
	print "Ready to do PSI analysis...\n";
	$starttime = time;
	my $adjpdefault = 0;
	$adjpdefault = 1 if($adjp == 1);
	my $commend = "perl " . $path . "/PSIsigma-PSI-v.1.1.pl " . $name . ".db " . $name . " " . $supporting_read_criteria . " " . $skipratio . " " . $intron_criteria . " " . $type . " " . $adjpdefault . " " . $denominator . " " . $ircheck . " " . $irrange;
	system($commend);
	$stoptime = time;
	$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
	$totaltime += $hours;
	print "===PSI analysis spent $hours hours.===\n";
	
	print "Filtering ΔPSI results...\n";
	$starttime = time;
	if($nfiles < 4 && $fmode == 0){
		print "Not enough samples for p-value calculation, so switch to fmode = 1.\n";
		$fmode = 1;
	}
	if($nfiles < 4 && $fmode == 2){
		print "Not enough samples for p-value calculation, so switch to fmode = 1.\n";
		$fmode = 1;
	}
	if($adjp == 2 && $fmode != 3){
		print "## Note: --fmode has been changed to 3 for p-value adjustment.\n";
		$fmode = 3;
	}
	$commend = "perl " . $path . "/PSIsigma-filter-v.1.0.pl " . $name . ".db " . $gtf . ".mapping.txt " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt " . $fmode;
	system($commend);
	$stoptime = time;
	$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
	$totaltime += $hours;
	print "===Filtering spent $hours hours.===\n";
	
	if(-e $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt"){
    	print "Archiving... " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt \n";
    	system("tar zcvf " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt.tar.gz " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt");
		system("rm " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".txt");
    }else{
    }
	
	if($adjp == 2){
		print "Adjusting p-value distribution...\n";
		$commend = "perl " . $path . "/PSIsigma-FDR-v.1.0.pl " . $name . "_r" . $supporting_read_criteria . "_ir" . $intron_criteria . ".sorted.txt 2 $trimp";
		system($commend);
		$stoptime = time;
		$hours = sprintf("%.4f",(($stoptime-$starttime)/3600));
		$totaltime += $hours;
		print "===P-value adjustment spent $hours hours.===\n";
	}

	
	print "\n***Total: $totaltime hours (or " . ($totaltime*60) . "mins).\n";
	
sub rundb{
	my $noveljunctioncriteria = shift;
	my $gtf = shift;
	my $chrs = shift;
	my $type = shift;
	my @chromosomes = split(/\t/,$chrs);
	foreach my $chromosome(@chromosomes){
		next if($chromosome=~/^chrGL/);
		next if($chromosome=~/^chrKI/);
		next if($chromosome=~/^GL/);
		next if($chromosome=~/^KI/);
		
		my $commend = "perl " . $path . "/PSIsigma-db-v.1.0.pl $gtf " . $chromosome . " " . $noveljunctioncriteria . " " . $type . " " . $irmode;
		#print "Doing... $commend\n";
		print "Doing... $chromosome\n";
		system("$commend");
	}
	
	system("cat *.db.tmp > $dbname");
	system("cat *.bed.tmp > $bedname");
	system("rm *.db.tmp");
	system("rm *.bed.tmp");
}

sub param{
	my $a = $_[0];
	my @array = @$a;
	
	if($array[0] eq "--help" || $array[0] eq "-h"){
		return "Help";
	}
	
	my %parameters;
	$parameters{"gtf"} = "-";
	$parameters{"name"} = "-";
	$parameters{"type"} = "-";
	$parameters{"nread"} = "-";
	$parameters{"fmode"} = "-";
	$parameters{"skipratio"} = "-";
	$parameters{"irmode"} = "-";
	$parameters{"adjp"} = "-";
	$parameters{"denominator"} = "-";
	$parameters{"irrange"} = "-";
	$parameters{"trimp"} = "-";
	
	my $oldformat = 1;
	for(my $i = 0;$i < scalar @array;$i++){
		if($array[$i]=~/^\-/){
			my $pam = $array[$i];
			$pam=~s/^\-\-//;
			$pam=~s/^\-//;
			if(!$parameters{$pam}){
				return "Can't recognize the parameter name: $pam";
			}else{
				if(!$array[($i+1)] && $array[($i+1)] != 0){
					return "Parameter $pam has no input value";
				}else{
					$oldformat = 0;
					$parameters{$pam} = $array[($i+1)];
				}
			}
		}
	}
	if($oldformat == 1){
		if($array[0]!~/\.gtf/){
			return "Not recognized parameter: $gtf";
		}
		return $array[0] . "\t" . $array[1] . "\t" . $array[2] . "\t" . $array[3] . "\t" . "0";
	}
	foreach my $key(keys %parameters){
		next if($key eq "fmode" || $key eq "skipratio" || $key eq "irmode" || $key eq "adjp" || $key eq "denominator" || $key eq "irrange" || $key eq "trimp");
		if($parameters{$key} eq "-"){
			print "Parameter $key has no input value";
			return "Parameter $key has no input value";
		}
	}
	if($parameters{"gtf"}!~/\.gtf/){
		return "--gtf parameter didn't find a files with .gtf extension";
	}
	if($parameters{"type"} != 1 && $parameters{"type"} != 2){
		return "--type parameter didn't find a correct number (1 or 2)";
	}
	return $parameters{"gtf"} . "\t" . $parameters{"name"} . "\t" . $parameters{"type"} . "\t" . $parameters{"nread"} . "\t" . $parameters{"fmode"} . "\t" . $parameters{"skipratio"} . "\t" . $parameters{"irmode"} . "\t" . $parameters{"adjp"} . "\t" . $parameters{"trimp"} . "\t" . $parameters{"denominator"} . "\t" . $parameters{"irrange"};
}

sub generateSJ{
	my $bam = shift;
	my $accession = shift;
	print "Generating .SJ.out.tab will need to re-sort $bam file by read names.\n";
	print "It will consume a lot of time, do you want to proceed? (Y/N)";
	#my $input = <STDIN>;
	my $input = "Y";
	chomp $input;
	my $newbam = "$accession.SortedbyName.bam";
	if($input ne "Y" && $input ne "N"){
		print "Bye.\n";
		exit;
	}else{
		print "\nYes, proceed.\n";
		my $accession = "-";
		my $dupcount = 0;
		open(INPUT, '-|',"samtools view " . $bam . " | head -n 100") or die $!;
		while (my $input = <INPUT>) {
			chomp $input;
			my @array = split(/\t/,$input);
			my ($name,$chr,$flag,$ss,$cigar) = ($array[0],$array[2],$array[1],$array[3],$array[5]);
			$accession = $name if($accession eq "-");
			$dupcount++ if($accession eq $name);
		}
		close(INPUT);
		if($dupcount > 90){
			print "[ERROR]: The .bam contains too many duplicated read names (over 90% in the first 100 lines).\n";
			exit;
		}
		print "Starting to sort $bam by read names...\n";
		#print "chromosome format = " . $chrformat . "\n";
		system("samtools sort -n $bam -o $newbam");
	}
	
	my $sjfn = $accession . ".SJ.out.tab";
	my $commend = "samtools view $newbam | awk -f " . $path . "/sjFromSAMcollapseUandM_inclOverlaps.awk > " . $sjfn;
	system($commend);
	#system("rm $newbam");
}