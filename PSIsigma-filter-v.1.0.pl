=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
	use strict;
    my ($db,$mapping,$input,$fmode) = @ARGV;

	print "Filtering mode = $fmode\n";
	#$fmode = 1;
	my $criteria = 0.75;
	my %symbol;
	my %strand;
	print "Reading... $mapping\n";
 	open(FILE,"$mapping") || die "Aborting.. Can't open $mapping : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        my ($ENST,$symbol,$strand) = split(/\t/,$line);
        $symbol{$ENST} = $symbol;
        $strand{$ENST} = $strand;
    }
    close(FILE);

	print "Reading... $db\n";
    my %ET;
    my %wing;
    my %names;
    open(FILE, "$db") || die "Aborting.. Can't open $db\n";
    while(my $line=<FILE>){
    	chomp $line;
    	next if($line eq "");
    	my ($chr,$i1s,$i1e,$i2s,$i2e,$tes,$tee,$anno,$as,$ae,$name,$gn) = split(/\t/,$line);
    	my ($et) = "-";
    	my ($accession,$num) = ($1,$2) if($name=~/(.*)\_(\d+)$/);

    	my @array = split(/\_/,$accession);
    	#$accession=~s/(.*)\_//;
    	my $ENST = $array[4];
    	if($ENST eq ""){
    		print "WARNNING: ENST format in the .db file is not correct\n";
    		print "accession = $accession\n";
    		$ENST = $accession;
    		$ENST=~s/(.*)\_//;
    		exit;
    	}
    	$wing{$name} = "$chr,$i1s,$i1e,$i2s,$i2e";
    	$ENST=~s/^Ex\.//;
    	$ENST=~s/^TSS\.//;
    	
    	if($name=~/\_W\_/){
    		if($i1s == $as && $i2e == $ae){
    			$ET{$accession} = "SES" if($num == 1);
    		}else{
    			$ET{$accession} = "MES" if($num > 1);
    		}
    	}
    	
    	if($name=~/\_S\_/){
    		if($i1s == $i1e){
    			$ET{$name} = "A3SS" if($strand{$ENST} eq "+");
    			$ET{$name} = "A5SS" if($strand{$ENST} eq "-");
    		}
    		if($i2s == $i2e){
    			$ET{$name} = "A5SS" if($strand{$ENST} eq "+");
    			$ET{$name} = "A3SS" if($strand{$ENST} eq "-");
    		}
    		if($i1s != $i1e && $i2s != $i2e){
    			print "[ERROR] coordinates are not A5SS or A3SS!\n";
    			print "$i1s,$i1e; $i2s,$i2e; $accession:$num\n";
    		}	
    	}
    	if($name=~/\_R\_/){
    		$ET{$name} = "IR";
    		$ET{$name} = "IR (overlapping region)" if($gn=~/\,/);
    	}
    }
    close(FILE);
    
	my %output;
	my %asae;
	my %es;
	my %ee;
	my %exonR;
	my ($maxn,$maxt) = (0,0);
 	open(FILE,"$input") || die "Aborting.. Can't open $input : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line=~/^Event ID/);
        my ($ID,$tmpgn,$exon,$eventtype,$N,$T,$exontype,$ENST,$dPSI,$pvalue,$FDR,$nvalues,$tvalues) = split(/\t/,$line);
        my @nvalue = split(/\|/,$nvalues);
        my @tvalue = split(/\|/,$tvalues);
        my $wing = $wing{$ID};
        if(!$wing{$ID}){
        	print "(ERROR) $ID can't find wings in the database.\n";
        	print "Aborting...\n";
        	exit;
        }
        if($fmode == 0){
        	next if($pvalue >= 0.01 || abs($dPSI) <= 10);
        }
        if($fmode == 1){
        	next if(abs($dPSI) <= 10);
        }
        if($fmode == 2){
        	next if($pvalue >= 0.05);
        }
        if($fmode == 3){
        	my ($avgn,$avgt) = (average(\@nvalue),average(\@tvalue));
        	#next if(abs($dPSI) < 5);
        	#next if($avgn < 10 && $avgt < 10);
        	#next if($avgn < 5 || $avgt < 5);
        	if($N > 1 && $T > 1){
        		#next if($pvalue > 0.05);
        	}
        	#my $no = 0;
        	#if($avgn < 20 || $avgt < 20){
        	#	next if($avgn < 20 && $avgt < 20);
        	#}else{
        	#	$no++;
        	#}
        	#if($avgn > 80 || $avgt > 80){
        	#	next if($avgn > 80 && $avgt > 80);
        	#}else{
        	#	$no++;
        	#}
			#if($no == 2){
			#	next if($pvalue >= 0.01 || abs($dPSI) <= 10);
			#}else{
        	#	next if($pvalue >= 0.5 || abs($dPSI) <= 5);
        	#}
        }
        $maxn = $N if($maxn < $N);
        $maxt = $T if($maxt < $T);
        $line=~s/\, /\|/g;
        my ($chr,$as,$ae,$t,$tmpENST,$num) = split(/\_/,$ID);
        $chr=~s/chr//;
        my ($exonchr,$exonss,$exonee) = split(/[\:\-]/,$exon);
        if($eventtype eq "R"){
        	$exon = "$chr:$as-$ae";
		}
		if(!$output{"$exon\t$wing\t$eventtype"}){
        	$output{"$exon\t$wing\t$eventtype"} = $line;
        }else{
        	my @tmparray = split(/\t/,$output{"$exon\t$wing\t$eventtype"});
        	my $tmpp = $tmparray[9];
        	next if($tmpp <= $pvalue);
        	$output{"$exon\t$wing\t$eventtype"} = $line;
        }
        if($eventtype eq "W"){
        	$es{$exonss}++;
        	$ee{$exonee}++;
        	$asae{"$chr\t$as\t$ae"}{"$exonss\t$exonee"}++;
        }
        if($eventtype eq "R"){
        	$exonR{$exon}++;
        }
    }
    close(FILE);
    
    $input=~s/\.txt$//;
    my $outfn = $input . ".sorted.txt";
    my $tmpfn = $input . ".filtered.txt";
    #$tmpfn = $input . ".0to1.filtered.txt" if($fmode == 2);
    #$outfn = $input . ".0to1.sorted.txt" if($fmode == 2);
    open(OUT,">" . $tmpfn) || die "Aborting.. Can't open " . $tmpfn . " : $!\n";
    print OUT "Event Region	Gene Symbol	Target Exon	Event Type	N	T	Exon Type	Reference Transcript	ΔPSI (%)	T-test p-value	FDR (BH)	N Values	T Values	Database ID\n";
    foreach my $exoninfo(keys %output){
    	my ($ID,$tmpgn,$exon,$eventtype,$N,$T,$exontype,$ENST,$dPSI,$pvalue,$FDR,$nvalues,$tvalues) = split(/\t/,$output{$exoninfo});
    	if($eventtype eq "R"){
    		if(!$exonR{$exon}){
    		}else{
    			#next if($exonR{$exon} > 1);
    		}
    	}
    	my ($chr,$as,$ae,$tmp,$tmpENST,$tmpnum) = split(/\_/,$ID);
    	#$chr=~s/chr//;
    	my ($accession,$num) = ($1,$2) if($ID=~/(.*)\_(\d+)$/);
    	#$accession=~s/(.*)\_//;

=begin		
    	if($ID=~/\_R1\./ || $ID=~/\_R2\./){
    		$ENST=~s/\./\_/g;
    		$ID=~s/\./\_/g;
    		($accession,$num) = ($1,$2) if($ID=~/(\w+)\_(\d+)$/);
    		#$accession=~s/\_/\./g;
    		$ENST=~s/\_/\./g;
    		#$accession=~s/(.*)\_R1/R1/;
    		#$accession=~s/(.*)\_R2/R2/;
    		#$accession=~s/\_/\./g;
    		#print "Name = $name\n";
    		#print "accession = $accession\n";
    		#print "num = $num\n";
    		#exit;
    	}
=end
=cut
    	if($fmode == 0){
    		next if($N < ($maxn*$criteria) || $T < ($maxt*$criteria));
    	}
    	my $symbol = "-";
    	$tmpENST=~s/^Ex\.//;
    	$tmpENST=~s/^TSS\.//;
    	if(!$symbol{$tmpENST}){
    		$symbol="N/A";
    		print "$tmpENST has no symbol!\n";
    	}else{
    		$symbol = $symbol{$tmpENST};
    	}
    	
    	my $ETID = $ID;
    	$ETID = $accession if($ID=~/\_W\_/);
    	
    	if(!$ET{$ETID}){
    		print "Can't find $accession\n";
    		print "ID = $ID\n";
    		print "ENST = $ENST\n";
    		print "symbol = $symbol\n";
    		exit;
    	}else{
    		$eventtype = $ET{$ETID};
    	}
    	if($eventtype eq "SES"){
    		my $mxs = 0;
    		my ($exonchr,$exonss,$exonee) = split(/[\:\-]/,$exon);
    		foreach my $loc(keys %{ $asae{"$chr\t$as\t$ae"} }){
    			my ($es,$ee) = split(/\t/,$loc);
    			next if($exonss <= $es && $exonee >= $es);
    			next if($exonss <= $ee && $exonee >= $ee);
    			next if($es < $exonss && $ee > $exonee);
    			$mxs = 1;
    			last;
    		}
    		$eventtype = "MXS" if($mxs == 1);
    	}
    	if($eventtype eq "A5SS" || $eventtype eq "A3SS"){
    		#my ($exonchr,$exonss,$exonee) = split(/[\:\-]/,$exon);
    		#if(!$es{$exonss} && !$ee{$exonee}){
    		#}else{
    		#	next if($fmode != 0);
    		#}
    	}
    	my $eventregion = "$chr:$as-$ae";
    	my ($exonchr,$exonss,$exonee) = split(/[\:\-]/,$exon);
    	if($eventtype=~/IR/){
    		$eventregion = $exon;
    		$exon = "$chr:$as-$ae";
    	}
    	if($ENST=~/^TSS/){
    		$eventtype = "TSS" . "|" . $eventtype;
    	}
    	#print OUT "$eventregion\t$symbol\t$exon\t$eventtype\t$N\t$T\t$exontype\t$ENST\t$dPSI\t$pvalue\t$FDR\t$nvalues\t$tvalues\n";
    	print OUT "$eventregion\t$symbol\t$exon\t$eventtype\t$N\t$T\t$exontype\t$ENST\t$dPSI\t$pvalue\t$FDR\t$nvalues\t$tvalues\t$ID\n";
    }
    close(OUT);
	
	my $sortcommand = "\(head \-n 1 " . $tmpfn . " \&\& tail \-n \+2 " . $tmpfn . " \| sort \-gr \-t '\t' \-k 9\) \> " . $outfn . "\n";
    system($sortcommand);
    system("rm " . $tmpfn);

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}