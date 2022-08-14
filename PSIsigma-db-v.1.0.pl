=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
	use strict;

    my ($gtf,$qchr,$noveljunctioncriteria,$longread,$irmode,$groupa,$groupb) = @ARGV;
    
    my ($gchr,$start,$end) = ("-",0,0);
    my $exons;
    my $count = 0;
    my $max = 0;
    
    my $starttime = time;

    my %anno;
    my %exonanno;
    open(FILE,"$gtf") || die "Aborting.. Can't open $gtf : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line=~/^\#/);
        my @array = split(/\t/,$line);
        next if($array[2] ne "gene" && $array[2] ne "exon");
        my ($chr,$cat,$start,$end,$strand,$name) = ($array[0],$array[1],$array[3],$array[4],$array[6],$array[8]);
        #$chr = "chr" . $chr if($chr!~/chr/);
        next if($chr ne $qchr);
        ##remove this for tomato
        #$name=~s/\.(\d+)\"\;/\"\;/g;
        if($array[2] eq "gene"){
        	$name=~s/(.*)gene\_name \"//;
        	$name=~s/\"\; (.*)//;
        	$name=~s/\"\;//;
        	$name=~s/gene\_id \"//;
        	if($name=~/\_/){
				$name=~s/\_/\./g;
			}
        	$anno{$chr}{"$start\t$end"}++;
        	$anno{"$chr"}{"$start\t$end"} = $name;
        }
        if($array[2] eq "exon"){
			my $ENST = $name;
 			$ENST=~s/(.*)transcript\_id \"//;
    		$ENST=~s/\"\; (.*)//;
    		$ENST=~s/\"\;//;
			if($ENST=~/\_/){
				$ENST=~s/\_/\./g;
			}
			
        	if($name=~/transcript_biotype/){
        		$name=~s/(.*)transcript_biotype \"//;
        		$name=~s/\"\; (.*)//;
        		$cat = $name;
        	}
        	if($name=~/transcript_type/){
				$name=~s/(.*)transcript_type \"//;
        		$name=~s/\"\; (.*)//;
        		$cat = $name;
        	}
        	my $label = -1;
        	$label = 1 if($cat eq "nonsense_mediated_decay");
        	#$exonanno{$chr}{"$start\t$end"}{$ENST}++;
        	if(!$exonanno{"$chr"}{"$start\t$end"}{$ENST}){
        		$exonanno{"$chr"}{"$start\t$end"}{$ENST} = $label;
        	}else{
        		$exonanno{"$chr"}{"$start\t$end"}{$ENST} = $label if($label == 1);
        	}
        }	
    }
    close(FILE);

    my @files;
 	open(FILE,"$groupa") || die "Aborting.. Can't open $groupa : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        my $accession = $line;
        $accession=~s/\.Aligned\.sortedByCoord\.out\.bam//;
		$accession=~s/\.sorted\.out\.bam//;
		$accession=~s/\.bam//;
		$accession=~s/\.$//;
		$accession.= ".SJ.out.tab";
        push(@files,$accession);
        
    }
    close(FILE);
 	open(FILE,"$groupb") || die "Aborting.. Can't open $groupb : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line eq "");
        my $accession = $line;
        $accession=~s/\.Aligned\.sortedByCoord\.out\.bam//;
		$accession=~s/\.sorted\.out\.bam//;
		$accession=~s/\.bam//;
		$accession=~s/\.$//;
		$accession.= ".SJ.out.tab";
        push(@files,$accession);
    }
    close(FILE);
  
    my %SJ;
    my %SJsite;
    my $countsj = 0;
    foreach my $jfn(@files){
    	next if(-z $jfn);
    	#print "$jfn...\n";
    	my $count = 0;
		open(FILE, "$jfn") || die "Aborting.. Can't open $jfn\n";
        while(my $line=<FILE>){
        	chomp $line;
        	my @array = split(/\t/,$line);
            my ($chr,$start,$end,$num) = ($array[0],$array[1],$array[2],$array[6]);
        	if(scalar @array == 7){
                #$jfn is a customized SJ file\n";
                $num = $array[3] if($longread == 1);
                $num = $array[4] if($longread == 2);
            }else{
                print "[ERROR:UNKNOWN FORMAT of $jfn]\n" if($longread == 2);
                $num = $array[7] if($longread == 2);
            }
            #$chr = "chr" . $chr if($chr!~/chr/);
            next if($chr ne $qchr);
            next if($num < $noveljunctioncriteria);
            $jfn=~s/\.SJ\.out\.tab//;
            $SJ{$chr}{"$start\t$end"}++;
            $SJsite{$chr}{$start}++;
            $SJsite{$chr}{$end}++;
            $count++;
        }
        close(FILE);
        $countsj++;
        #print "number of valid junction = $count\n";
    }
    print "Number of valid .SJ.out.tab files = $countsj\n";

    my %tmp;
    my %bed;
    my %db;
    my %unique;
    my $total_exon = 0;
    open(FILE,"$gtf") || die "Aborting.. Can't open $gtf : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        my $num = scalar keys %tmp;
        #if(eof){
        #}else{
        #	next if($line!~/^$qchr/);
        #}
        next if($line!~/\texon\t/ && $line!~/\ttranscript\t/);
        my @array = split(/\t/,$line);
        my ($chr,$cat,$ss,$ee,$strand,$anno) = ($array[0],$array[2],$array[3],$array[4],$array[6],$array[8]);
		#$chr = "chr" . $chr if($chr!~/chr/);
		next if($chr ne $qchr);
        $gchr = $chr if($gchr eq "-");
        if($start == 0 && $end == 0 && $cat eq "transcript"){
        	$gchr = $chr;
        	$start = $ss;
        	$end = $ee;
        }

		if($cat eq "exon" && $gchr eq $chr){
			my $tid = "";
			my $tid = $anno;
			$tid=~s/(.*)transcript\_id \"//;
    		$tid=~s/\"\; (.*)//;
    		$tid=~s/\"\;//;
			if($tid=~/\_/){
				$tid=~s/\_/\./g;
			}
       		push(@{$tmp{$tid}},"$chr\t$ss\t$ee");
       		if(!$unique{"$chr\t$ss\t$ee"}){
       			$unique{"$chr\t$ss\t$ee"}++;
       			$total_exon++;
       		}
		}
		
		#strand is not needed because exons were inserted to @tmp array based on their coordinates instead of exon number
        if($cat eq "transcript" || eof){
        	if($ss > $end || eof){
        		$count++;
        		my $num = scalar keys %tmp;
				if($num >= 1){
					$max = $num if($num > $max);
					my $id = $gchr . "_" . $start . "_" . $end . "_" . $strand;
					run($id, $total_exon, %tmp);
				}
        		%tmp = ();
        		%unique = ();
        		$gchr = $chr;
        		$start = $ss;
        		$end = $ee;
        		$total_exon = 0;
        	}else{
        		if($ee > $end && $gchr eq $chr){
        			$end = $ee;
				}
			}
			next;
		}
        

	}

    my $stoptime = time;
    my $mins = sprintf("%.2f",(($stoptime-$starttime) / 60));
    
	open(OUT,">$qchr.bed.tmp") || die "Aborting.. Can't open $qchr.bed.tmp : $!\n";
	foreach my $bedoutput(sort keys %bed){
		print OUT $bed{$bedoutput} . "\n";
	}
	close(OUT);
	
	open(OUT,">$qchr.db.tmp") || die "Aborting.. Can't open$qchr.db.tmp : $!\n";
	foreach my $dboutput(sort keys %db){
		print OUT $db{$dboutput} . "\n";
	}
	close(OUT);
    
    print "time = $mins mins\n";
	open(OUT,">>time.txt") || die "Aborting.. Can't open time.txt : $!\n";
	print OUT "$qchr\t" . ($mins) . " mins\n";
	close(OUT);

sub run2
{
	print "run2\n";
	return 0;
}
sub run
{
	my ($id,$total_exon,%input) = @_;
	my ($ichr,$istart,$iend,$istrand) = split(/\_/,$id);

	my $starttimex = time;
	my $testcount = scalar keys %input;

	if($testcount > 10000){
		my $showcount = 0;
		foreach my $tid(sort keys %input){
			print "$tid\n";
			$showcount++;
			last if($showcount >= 5);
		}
	}
	
	#my $fn = "count_exon.txt";
	#if($total_exon > 500){
	#	open(OUT,">>$fn") || die "Aborting.. Can't open $fn : $!\n";
	#	print OUT "$id\t$total_exon\n";
	#	close(OUT);
	#}

	my %intron;
	
	my %knownss;
	my %knownes;
	my %knownexons;
		
	foreach my $tid(sort keys %input){
		next if($tid eq "");
		my @array = @{$input{$tid}};
		my $num = scalar @array;
		for(my $i = 1;$i < $num;$i++){
			my ($chr1,$ss1,$ee1) = split(/\t/,$array[$i-1]);
			my ($chr2,$ss2,$ee2) = split(/\t/,$array[$i]);
			$knownexons{"$ss1\t$ee1"}++;
			$knownexons{"$ss2\t$ee2"}++;
			$knownss{($ee1+1)}{"$ss2,$ee2"}++;
			$knownes{($ss2-1)}{"$ss1,$ee1"}++;
			$intron{($ee1+1) . "\t" . ($ss2-1)}++;
		}
	}
	
	#Add mutually exclusive introns
	my @knownss;
	foreach my $site(sort keys %knownss){
		next if(scalar keys %{$knownss{$site}} == 1);
		push(@knownss,$site);
	}
	my @knownes;
	foreach my $site(sort keys %knownes){
		next if(scalar keys %{$knownes{$site}} == 1);
		push(@knownes,$site);
	}
	foreach my $ss(@knownss){
		foreach my $es(@knownes){
			my $bind = 0;
			foreach my $bindexon1(keys %{$knownss{$ss}}){
				foreach my $bindexon2(keys %{$knownes{$es}}){
					if($bindexon1 eq $bindexon2){
						$bind++;
					}
				}
				last if($bind > 1);
			}
			$intron{$ss . "\t" . $es}++ if($bind > 1);
		}
	}

	## Seek for novel exons
	my %novelintron;
	my %SJsites;
	foreach my $sjintron(sort keys %{ $SJ{$ichr} }){
		my ($is,$ie) = split(/\t/,$sjintron);
		next if($is < $istart || $ie > $iend);
		
		if(!$intron{$sjintron}){
			$novelintron{$sjintron}++;
			$intron{$sjintron}++;
		}else{
			next;
		}
	}

	my %collection;
	my %queue;
	my %last;
	foreach my $intron(sort keys %intron){
		my ($is,$ie) = split(/\t/,$intron);
		my ($targetexon,$exonanno,$wings);
		foreach my $tid(sort keys %input){
			next if($tid eq "");
			my @array = @{$input{$tid}};
			my $num = scalar @array;
			my ($lastchr,$lastss,$lastee) = split(/\t/,$array[$num-1]);
			if(($lastss) < $is && ($lastee) > $ie){
				$targetexon = ($lastss) . "\t" . ($lastee);
				$exonanno = $exonanno{"$lastchr"}{"$lastss\t$lastee"}{$tid};
				$wings = ($is-$is) . "\t" . ($is-$is) . "\t" . ($is-$is) . "\t" . ($ie-$is);
				$targetexon = ($lastss-$is) . "\t" . ($lastee-$is);
				if(!$queue{"$lastchr\t$is\t$ie\tR\t" . $wings . "\t" . $targetexon}){
					$queue{"$lastchr\t$is\t$ie\tR\t" . $wings . "\t" . $targetexon}++;
				}else{
					next;
				}
				if($irmode == 1){
					$collection{"$lastchr\t$is\t$ie\tR\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
				}else{
					if(!$novelintron{$intron}){
						$collection{"$lastchr\t$is\t$ie\tR\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					}
				}
				last;
			}
			next if($lastee < $is);
			for(my $i = 0;$i < ($num-1);$i++){
				my ($chr2,$ss2,$ee2) = split(/\t/,$array[$i]);
				my ($chr3,$ss3,$ee3) = split(/\t/,$array[$i+1]);
				
				last if($is == ($ee2+1) && $ie == ($ss3-1));
				next if($is >= ($ee3+1));
				last if($ie <= ($ss2-1));
				
				#For regular intron retention event.
				if(($ss2-1) < $is && $ie < ($ee2+1)){
					$targetexon = ($ss2) . "\t" . ($ee2);
					$exonanno = $exonanno{$chr2}{$targetexon}{$tid};
					$targetexon = ($ss2-$is) . "\t" . ($ee2-$is);
					$wings = ($is-$is) . "\t" . ($is-$is) . "\t" . ($is-$is) . "\t" . ($ie-$is);
					if(!$queue{"$lastchr\t$is\t$ie\tR\t" . $wings . "\t" . $targetexon}){
						$queue{"$lastchr\t$is\t$ie\tR\t" . $wings . "\t" . $targetexon}++;
					}else{
						next;
					}
					if($irmode == 1){
						$collection{"$ichr\t$is\t$ie\tR\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					}else{
						if(!$novelintron{$intron}){
							$collection{"$ichr\t$is\t$ie\tR\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
						}
					}
					last;
				}
				
                #For novel IR events (aggressive mode)
                if($irmode == 1){
                	my ($newis,$newie) = (($ee2+1),($ss3-1));
                	my $newtid = "Ex." . $tid;
                	my $newtargetexon = ($ss2-$newis) . "\t" . ($ee3-$newis);
                	$wings = ($newis-$newis) . "\t" . ($newis-$newis) . "\t" . ($newis-$newis) . "\t" . ($newie-$newis);
                	$exonanno = 2;
                	my $knownIR = 0;
                	foreach my $ckey(keys %collection){
                    	my ($cchr,$cis,$cie,$ctype,$cid) = split(/\t/,$ckey);
                    	my ($x1sf,$x1ef,$x2sf,$x2ef,$xtess1,$xteee1,$xeanno1) = split(/\t/,$collection{$ckey});
                    	next if($ctype ne "R");
                    	if($xtess1 == ($ss2-$newis) && $xteee1 == ($ee2-$newis)){
                        	$knownIR = 1;
                        	last;
                    	}
                    	if(($cis == $newis) && ($cie == $newie)){
                        	$knownIR = 1;
                        	last;
                    	}
                	}
                	if($knownIR == 0){
                		if(!$queue{"$ichr\t$newis\t$newie\tR\t" . $wings . "\t" . $newtargetexon}){
							$queue{"$ichr\t$newis\t$newie\tR\t" . $wings . "\t" . $newtargetexon}++;
						}else{
							next;
						}
                    	$collection{"$ichr\t$newis\t$newie\tR\t$newtid"} .= "," . $wings . "\t" . $newtargetexon . "\t" . $exonanno;
                	}
				}
                last if($is == ($ee2+1) && $ie == ($ss3-1));
				
				#For novel ASS events
				if(!$novelintron{$intron}){
				}else{
					#left site the same, new splice site in the intron
					if(($ee2+1) == $is && $ie < ($ss3-1)){
						my $tmpintron = ($is) . "\t" . ($ss3-1);
						$targetexon = ($ie+1-$is) . "\t" . ($ss3-$is);
						$exonanno = 2;
						$wings = ($is-$is) . "\t" . ($is-$is) . "\t" . ($is-$is) . "\t" . ($ie-$is);
						my $TSS = 0;
						if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
							$TSS++;
						}
						if($i == 0 || $i == ($num-2)){
							$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
						}
						if(!$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$tid = "Ex." . $tid;
						$collection{"$ichr\t$tmpintron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
						last;
					}
					#left site the same, new splice site in the exon
					if(($ee2+1) == $is && $ie > ($ss3-1) && $ie < $ee3){
						my $tmpintron = ($is) . "\t" . ($ie);
						$targetexon = ($ss3-$is) . "\t" . ($ie-$is);
						$exonanno = 2;
						$wings = ($is-$is) . "\t" . ($is-$is) . "\t" . ($is-$is) . "\t" . (($ss3-1)-$is);
						my $TSS = 0;
						if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
							$TSS++;
						}
						if($i == 0 || $i == ($num-2)){
							$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
						}
						if(!$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$tid = "Ex." . $tid;
						$collection{"$ichr\t$tmpintron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
						last;
					}
					#right site the same, new splice site in the intron
					if(($ee2+1) < $is && $ie == ($ss3-1)){
						my $tmpintron = ($ee2+1) . "\t" . ($ie);
						$targetexon = (($ee2+1)-($ee2+1)) . "\t" . ($is-($ee2+1)-1);
						$exonanno = 2;
						$wings = ($is-($ee2+1)) . "\t" . ($ie-($ee2+1)) . "\t" . ($ie-($ee2+1)) . "\t" . ($ie-($ee2+1));
						my $TSS = 0;
						if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
							$TSS++;
						}
						if($i == 0 || $i == ($num-2)){
							$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
						}
						if(!$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$tid = "Ex." . $tid;
						$collection{"$ichr\t$tmpintron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
						last;
					}
					#right site the same, new splice site in the exon
					if(($ee2+1) > $is && $ie == ($ss3-1) && $ss2 < $is){
						my $tmpintron = ($is) . "\t" . ($ie);
						$targetexon = ($is-($is)) . "\t" . ($ee2-($is));
						$exonanno = 2;
						$wings = (($ee2+1)-($is)) . "\t" . ($ie-($is)) . "\t" . ($ie-($is)) . "\t" . ($ie-($is));
						my $TSS = 0;
						if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
							$TSS++;
						}
						if($i == 0 || $i == ($num-2)){
							$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
						}
						if(!$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$tmpintron\tS\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$tid = "Ex." . $tid;
						$collection{"$ichr\t$tmpintron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
						last;
					}
				}

				#For known alternative splice site event.
				if(($ee2+1) == $is && $ie > ($ss3-1) && $ie < ($ee3+1)){
					#$targetexon = ($ss3) . "\t" . ($ie);
					$exonanno = $exonanno{$chr3}{$targetexon}{$tid};
					$targetexon = ($ss3-$is) . "\t" . ($ie-$is);
					$wings = ($is-$is) . "\t" . ($is-$is) . "\t" . ($is-$is) . "\t" . (($ss3-1)-$is);
					my $TSS = 0;
					if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
						$TSS++;
					}
					if($i == 0 || $i == ($num-2)){
						$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
					}
					if(!$queue{"$ichr\t$intron\tS\t" . $wings . "\t" . $targetexon}){
						$queue{"$ichr\t$intron\tS\t" . $wings . "\t" . $targetexon}++;
					}else{
						next;
					}
					$collection{"$ichr\t$intron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					last;
				}
				if(($ee2+1) > $is && $ie == ($ss3-1) && $is > ($ss2-1)){
					#$targetexon = ($is) . "\t" . ($ee2);
					$exonanno = $exonanno{$chr2}{$targetexon}{$tid};
					$targetexon = ($is-$is) . "\t" . ($ee2-$is);
					$wings = (($ee2+1)-$is) . "\t" . ($ie-$is) . "\t" . ($ie-$is) . "\t" . ($ie-$is);
					my $TSS = 0;
					if(!$SJsite{$chr2}{($ss2-1)} || !$SJsite{$chr2}{($ee3+1)}){
						$TSS++;
					}
					if($i == 0 || $i == ($num-2)){
						$tid = "TSS." . $tid if($num > 2 && $TSS == 1);
					}
					if(!$queue{"$ichr\t$intron\tS\t" . $wings . "\t" . $targetexon}){
						$queue{"$ichr\t$intron\tS\t" . $wings . "\t" . $targetexon}++;
					}else{
						next;
					}
					$collection{"$ichr\t$intron\tS\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					last;
				}
				
				#exon skipping for 1 or multiple
				next if($i == 0);
				my ($chr1,$ss1,$ee1) = split(/\t/,$array[$i-1]);
				last if($is == ($ee1+1) && $ie == ($ss2-1));
				if($is <= ($ee1+1) && $ie >= ($ss3-1)){
					if(!$collection{"$ichr\t$intron\tW\t$tid"}){
						next if(($ee1+1) != $is);
						$targetexon = ($ss2) . "\t" . ($ee2);
						$exonanno = $exonanno{$chr2}{$targetexon}{$tid};
						$targetexon = ($ss2-$is) . "\t" . ($ee2-$is);
						$wings = ($ee1+1-$is) . "\t" . ($ss2-1-$is) . "\t" . ($ee2+1-$is) . "\t" . ($ss3-1-$is);
						if(!$exonanno){
							print "$tid [$chr2,$ss2,$ee2] can't find exon type.\n";
							exit;
						}
						if(!$queue{"$ichr\t$intron\tW\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$intron\tW\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$collection{"$ichr\t$intron\tW\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					}else{
						$targetexon = ($ss2) . "\t" . ($ee2);
						$exonanno = $exonanno{"$chr2"}{$targetexon}{$tid};
						$targetexon = ($ss2-$is) . "\t" . ($ee2-$is);
						$wings = ($ee1+1-$is) . "\t" . ($ss2-1-$is) . "\t" . ($ee2+1-$is) . "\t" . ($ss3-1-$is);
						if(!$exonanno){
							print "$tid [$chr2,$ss2,$ee2] can't find exon type.\n";
							exit;
						}
						if(!$queue{"$ichr\t$intron\tW\t" . $wings . "\t" . $targetexon}){
							$queue{"$ichr\t$intron\tW\t" . $wings . "\t" . $targetexon}++;
						}else{
							next;
						}
						$collection{"$ichr\t$intron\tW\t$tid"} .= "," . $wings . "\t" . $targetexon . "\t" . $exonanno;
					}
				}

			}
		}
	}	
	
		#Two S events can create a novel exon if they have the same constitutive exons
		#Doesn't work. block now!
=begin
		my %novelexon;
		my $novel = 0;
		my $create = 0;
		foreach my $event1(keys %collection){
			next if($event1!~/\tS\t/);
			my $junctions = $collection{$event1};
			$junctions=~s/\,//;
			my @junctions = split(/\,/,$junctions);
			next if(scalar @junctions == 1);
			my $num = scalar @junctions;
			#Left
			my ($js11,$je11,$js12,$je12) = split(/\t/,$junctions[0]);
			$js12 = $js11 if($js12 == $je12);
			$je12 = $je11 if($js12 == $je12);
			#Right
			my ($js21,$je21,$js22,$je22) = split(/\t/,$junctions[($num-1)]);
			
			next if($js12 == $js22 || $je12 == $je22);
			my ($chr1,$is1,$ie1,$type1,$tid1) = split(/\t/,$event1);
			$tid1=~s/^TSS\.//;
			$novelexon{$tid1}{($je12+1) . "\t" . ($js22-1)}++;
			undef $collection{$event1};
		}
			
	foreach my $tid(keys %novelexon){
		#print "checking... novel exon for $tid\n";
		my ($targetexon,$exonanno);
		foreach my $loc(keys %{ $novelexon{$tid} }){
			my ($tes,$tee) = split(/\t/,$loc);
			if(!$input{$tid}){
				$tid=~s/^TSS\.//;
			}
			my @array = @{$input{$tid}};
			my $num = scalar @array;
			for(my $i = 0;$i < ($num-1);$i++){
				my ($chr2,$ss2,$ee2) = split(/\t/,$array[$i]);
				next if($ss2 > $tes);
				my ($chr3,$ss3,$ee3) = split(/\t/,$array[$i+1]);
				if($ee2 < $tes && $ss3 > $tee){
					my $intron = ($ee2+1) . "\t" . ($ss3-1);
					print "Possible false novel exon skipping event " . ($ee2+1) . "-" . ($tes-1) . "\t" . ($tee+1) . "-" . ($ss3-1) . "\n" if(($ee2+1) == ($tes-1));
					$targetexon = ($tes-($ee2+1)) . "\t" . ($tee-($ee2+1));
					$exonanno = 2;
					$collection{"$ichr\t$intron\tW\t$tid"} .= "," . ($ee2+1-($ee2+1)) . "\t" . ($tes-1-($ee2+1)) . "\t" . ($tee+1-($ee2+1)) . "\t" . ($ss3-1-($ee2+1)) . "\t" . $targetexon . "\t" . $exonanno;
					last;
				}
			}
		}
	}
=end
=cut
	
	#Remove duplicates
	#Don't need this anymore. Block now.
=begin
	my %empty;
	foreach my $a1(sort keys %collection){
		next if(!$collection{$a1});
		if(!$empty{$a1}){
		}else{
			next;
		}
		my ($achr,$as,$ae,$acat,$atid) = split(/\t/,$a1);
		my @a1 = split(/\t/,$collection{$a1});
		my $labela1 = $a1[0] . "\t" . $a1[1] . "\t" . $a1[2] . "\t" . $a1[3];
		foreach my $a2(sort keys %collection){
			next if($a1 eq $a2);
			next if(!$collection{$a2});
			my ($bchr,$bs,$be,$bcat,$btid) = split(/\t/,$a2);
			next if($acat ne $bcat);
			my @a2 = split(/\t/,$collection{$a2});
			my $labela2 = $a2[0] . "\t" . $a2[1] . "\t" . $a2[2] . "\t" . $a2[3];
			if($collection{$a1} eq $collection{$a2} && $acat ne "S"){
				$empty{$a2} = 1;
				next;
			}
			if($labela1 eq $labela2 && $acat eq "S"){
				$empty{$a2} = 1;
				next;
			}
			
		}
	}
	
	foreach my $accession(sort keys %collection){
		my ($achr,$as,$ae,$acat,$atid) = split(/\t/,$accession);
		if(!$empty{$accession}){
		}else{
			undef $collection{$accession};
			next;
		}
	}
=end
=cut

    my $bedoutput;
    my $dboutput;
    foreach my $c(sort keys %collection){
        next if(!$collection{$c});
        my @groups = split(/\,/,$collection{$c});
        delete $collection{$c};
        my $last = scalar @groups;
        my ($ichr,$iss,$iee,$icat,$tid) = split(/\t/,$c);
        for(my $i = 1;$i < $last;$i++){
        	my @tmp = split(/\t/,$groups[$i]);
        	my $newvalue;
        	my $numtmp = scalar @tmp;
			for(my $j = 0;$j < $numtmp;$j++){
				###change here for tomato
				#if($j < ($numtmp-1)){
				if($j < ($numtmp)){
					if(!$newvalue){
						$newvalue = ($tmp[$j]+$iss);
					}else{
						$newvalue .= "\t" . ($tmp[$j]+$iss);
					}
				}else{
					$newvalue .= "\t" . ($tmp[$j]);
				}
			}
			$groups[$i] = $newvalue;
		}
        my $accession = $c;
        $accession=~s/\t/\_/g;
        my $names = "-";
        foreach my $loc(sort keys %{ $anno{$ichr} }){
        	my ($start,$end) = split(/\t/,$loc);
        	next if($end < $iss || $start > $iss);
        	$names .= ", " . $anno{$ichr}{$loc};
        }
        $names=~s/\-\, //;
        
        
		###Check W connectivity
        if($icat eq "W"){
        	my $disconnect = 0;
        	for(my $k = 1;$k < ($last-1);$k++){
        		next if(!$groups[$k]);
        		my ($intron1sf,$intron1ef,$intron2sf,$intron2ef,$tess1,$teee1,$eanno1) = split(/\t/,$groups[$k]);
        		my ($intron1sl,$intron1el,$intron2sl,$intron2el,$tess2,$teee2,$eanno2) = split(/\t/,$groups[$k+1]);
        		$disconnect = 1 if($intron2sf != $intron1sl || $intron2ef != $intron1el);
        		last if($intron2sf != $intron1sl || $intron2ef != $intron1el);
			}
			my ($intron1sf,$intron1ef,$intron2sf,$intron2ef,$tess1,$teee1,$eanno1) = split(/\t/,$groups[1]);
        	my ($intron1sl,$intron1el,$intron2sl,$intron2el,$tess2,$teee2,$eanno2) = split(/\t/,$groups[$last-1]);

        	$disconnect = 1 if($iss != $intron1sf || $iee != $intron2el);
			next if($disconnect == 1);
    	}
    	
    	#####DB output#####################################
        for(my $i = 1; $i < scalar @groups; $i++){
        	my ($intron1sf,$intron1ef,$intron2sf,$intron2ef,$tess,$teee,$eanno) = split(/\t/,$groups[$i]);
        	$eanno = $exonanno{$ichr}{"$tess\t$teee"}{$tid};
        	$eanno = "-" if($eanno eq "-1");
        	$eanno = "NMD" if($eanno eq "1");
        	$eanno = "Novel" if($eanno eq "2");
        	if($exonanno{$ichr}{"$tess\t$teee"}{$tid} eq ""){
        		#print "$accession [$ichr:$tess-$teee] of $tid has no exon type\n";
        		$eanno = "-";
        	}
        	$dboutput .= "$ichr\t$intron1sf\t$intron1ef\t$intron2sf\t$intron2ef\t$tess\t$teee\t$eanno\t$iss\t$iee\t$accession" . "_" . "$i\t$names\n";
		}
        
        ####BED output######################################
        
        my ($a,$b) = ("","");
        for(my $i = 1; $i < scalar @groups; $i++){
        if($icat eq "S"){
        	$bedoutput .= "\r$ichr\t$iss\t$iee\t$accession\t0\t+\t$iss\t$iee\t255,0,0\t" . ($last);
			my @tmp = split(/\t/,$groups[$i]);
			if($tmp[0] == $tmp[1]){
				$a = ($tmp[3]-$tmp[2]+1) . "," . "1";
				$b = "0" . "," . ($iee-$iss+1);
			}else{
				$a = "1," . ($tmp[1]-$tmp[0]+1);
				$b = "0," . ($tmp[0]-$iss+1);
			}
		}
		if($icat eq "R"){
			my @tmp = split(/\t/,$groups[$i]);
			$bedoutput .= "\r$ichr\t$iss\t" . $tmp[3] . "\t$accession\t0\t+\t$iss\t" . $tmp[3] . "\t255,0,0\t2";
			if($tmp[2] < $iss){
				$a = ($tmp[3]-$tmp[2]+1) . ",0";
				$b = "0," . ($iss-$tmp[2]);
			}else{
				$a = "0," . ($tmp[3]-$tmp[2]+1);
				$b = "0," . ($tmp[2]-$iss);
			}
		}
		}
		
		if($icat eq "W"){
			$b = "0";
			$bedoutput .= "\r$ichr\t$iss\t$iee\t$accession\t0\t+\t$iss\t$iee\t255,0,0\t" . ($last);
			for(my $i = 1;$i < $last;$i++){
				my @tmp = split(/\t/,$groups[$i]);
				$a .= "," . ($tmp[1]-$tmp[0]+1);
				$a .= "," . ($tmp[3]-$tmp[2]+1) if($i == ($last-1));
				$b .= "," . ($tmp[2]-$iss+1) if($i < $last);
			}
			$a=~s/\,//;
		}
		$bedoutput .= "\t" . $a . "\t" . $b;
		
	}
	my $cc = scalar keys %collection;
	undef %input;
	undef %collection;
	if(!$bedoutput){
	}else{
		$bed{$id} = $bedoutput;
		$db{$id} = $dboutput;
	}

	my $stoptimex = time;
    my $seconds = sprintf("%.4f",($stoptimex-$starttimex));
    
	return $seconds;
}



#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################


