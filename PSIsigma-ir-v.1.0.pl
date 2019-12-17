=begin
PSI-Sigma: A splicing-detection method for short-read and long-read RNA-seq data
© Kuan-Ting Lin, 2018-2024
PSI-Sigma is free for non-commercial purposes by individuals at an academic or non-profit institution.
For commercial purposes, please contact tech transfer office of CSHL via narayan@cshl.edu
=end
=cut
#!/usr/bin/perl -w
	use strict;
 	use Statistics::Basic qw(:all);
	
    my ($db,$bam,$type) = @ARGV;
    
    $bam=~s/(.*)\///;
    
    my $genome = "Human";
    if($bam=~/Sequins/){
    	$genome = "IS";
    }
    my %tmp;
    my %unique;
    my ($gchr,$start,$end) = ("-",0,0);
    my $exons;
    my $max = 0;
    my %bed;
    
    my $starttime = time;

    my %bin;
    my $binsize = 100;
    my @chromosomes;
    my $chrformat = 0;
    my $readformat = "single-end";
    
	if($genome eq "IS"){
		push(@chromosomes,"IS");
		$chrformat = 1;
	}else{
		for(my $i = 1; $i <= 22; $i++){
			push(@chromosomes,$i);
		}
		push(@chromosomes,"X");
		push(@chromosomes,"Y");
		push(@chromosomes,"M");
		push(@chromosomes,"MT");
		open(INPUT, '-|',"samtools view " . $bam . " | head -n 1") or die $!;
		while (my $input = <INPUT>) {
			chomp $input;
			my @array = split(/\t/,$input);
			my ($name,$chr,$flag,$ss,$cigar) = ($array[0],$array[2],$array[1],$array[3],$array[5]);
			if($chr=~/^chr/){
				$chrformat = 1;
			}
		}
		close(INPUT);
		print "chromosome format = " . $chrformat . "\n";
	}
	
	open(INPUT, "-|","samtools view -f 2 " . $bam . " | head -n 5") or die $!;
	while (my $input = <INPUT>) {
		chomp $input;
		$readformat = "paired-end";
	}
	close(INPUT);
	print "read format = " . $readformat . "\n";
	if($readformat eq "single-end"){
		print "[Warnning]: The use of single-end short reads is not recommended for alternative splicing analysis.\n"; 
	}
	
	system("mkdir \_wososatmp");

foreach my $ichr(@chromosomes){
	my %anno;
    my %introns;
    print "Doing... chromosome $ichr\n";
    open(FILE,"$db") || die "Aborting.. Can't open $db : $!\n";
    while(my $line=<FILE>){
        chomp $line;
        next if($line=~/^\#/);
        next if($line eq "");
		my ($chr,$i1s,$i1e,$i2s,$i2e,$tes,$tee,$anno,$as,$ae,$name,$gn) = split(/\t/,$line);
        $chr=~s/chr//;
		next if($chr ne $ichr);
        next if($name!~/\_R\_/);
		my $c1 = $as + 1;
		my $c2 = $ae - 1;
		my $interval = sprintf("%d",($ae-$as)/5);
		my $c3 = ($as+$interval-1);
		my $c4 = ($as+($interval*2)-1);
		my $c5 = ($as+($interval*3)-1);
		my $loc = "$i2s\t$i2e";
        $anno{$chr}{$c1}{$loc}++;
        $anno{$chr}{$c2}{$loc}++;
        $anno{$chr}{$c3}{$loc}++;
        $anno{$chr}{$c4}{$loc}++;
        $anno{$chr}{$c5}{$loc}++;   
        $introns{$chr}{$loc}{$c1} = 0;
        $introns{$chr}{$loc}{$c2} = 0;
        $introns{$chr}{$loc}{$c3} = 0;
        $introns{$chr}{$loc}{$c4} = 0;
        $introns{$chr}{$loc}{$c5} = 0;
    }
    close(FILE);

	my $count = 0;
	my $add = "";
	##For Illumina short reads
	if($type == 1){
		$add = "-F 256 -f 2 -q 255" if($readformat eq "paired-end");
		$add = "-F 256 -q 255" if($readformat eq "single-end");
	}
	if($chrformat == 1 && $ichr!~/^chr/){
		$ichr = "chr" . $ichr;
	}
	if($chrformat == 0){
		$ichr=~s/chr//;
	}
	open(INPUT, '-|',"samtools view " . $add . " " . $bam . " " . $ichr . "") or die $! . "\r" . "(error)$add $bam $ichr\r";
	while (my $input = <INPUT>) {
		chomp $input;
		next if($input eq "");
		$count++;
		my @array = split(/\t/,$input);
		my ($name,$chr,$flag,$ss,$cigar) = ($array[0],$array[2],$array[1],$array[3],$array[5]);
		if(!$cigar){
			print "problematic input = $input\n";
			next;
		}
		##For Illumina short reads
		if($type == 1){
			next if($cigar=~/[NSH*]/);
		}
		##For MinION long reads
		if($type == 2){
			next if($cigar eq "*");
		}
		$chr=~s/chr//;
		my @C = split(/[A-Z]/,$cigar);
		my @L = split(/[0-9]*/,$cigar);
		my $nc = scalar @C;
		my $nl = scalar @L;
		
		if($type == 1){
			for(my $i=0;$i<$nc;$i++){
				if(!$L[$i+1]){
					print "i = $i\n";
					print "c = $nc\n";
					print "l = $nl\n";
					print "C[$i] = " . $C[$i] . "\n";
					print "L[$i] = " . $L[$i] . "\n";
					print "cigar = $cigar\n";
					exit;
				}
				if($L[$i+1] ne "M"){
					$ss+=$C[$i] if($L[$i+1] eq "N" || $L[$i+1] eq "D");
				}else{
					my $nb = $C[$i];
					for(my $j=0;$j<$nb;$j++){
						my $pos = $ss+$j;
						if(scalar keys %{$anno{$chr}{$pos}} == 0){
							next;
						}else{
							foreach my $loc(keys %{ $anno{$chr}{$pos} }){
								next if(!$anno{$chr}{$pos}{$loc});
								$introns{$chr}{$loc}{$pos}++;
							}
						}
					}
					$ss+=$C[$i];
				}
			}
		}
		if($type == 2){
			my $longN = 0;
			for(my $i=0;$i<$nc;$i++){
				if($L[$i+1] eq "N"){
					if($C[$i] > 50){
						$longN = 1;
						last;
					}
				}
			}
			next if($longN == 1);
			for(my $i=0;$i<$nc;$i++){
				if(!$L[$i+1]){
					print "i = $i\n";
					print "c = $nc\n";
					print "l = $nl\n";
					print "C[$i] = " . $C[$i] . "\n";
					print "L[$i] = " . $L[$i] . "\n";
					print "cigar = $cigar\n";
					exit;
				}
				if($L[$i+1] eq "N" || $L[$i+1] eq "I" || $L[$i+1] eq "S"){
					$ss+=$C[$i] if($L[$i+1] eq "N");
				}else{
					my $nb = $C[$i];
					for(my $j=0;$j<$nb;$j++){
						my $pos = $ss+$j;
						if(scalar keys %{$anno{$chr}{$pos}} == 0){
							next;
						}else{
							foreach my $loc(keys %{ $anno{$chr}{$pos} }){
								next if(!$anno{$chr}{$pos}{$loc});
								$introns{$chr}{$loc}{$pos}++;
							}
						}
					}
					$ss+=$C[$i];
				}
			}
		}
	}
	close(INPUT);

    if($ichr!~/^chr/){
    	$ichr = "chr" . $ichr;
    }
    
	open(OUT,">_wososatmp/" . $bam . "." . $ichr . ".txt") || die "Aborting.. Can't open outputtmpfile.txt : $!\n";
	foreach my $chr(sort keys %introns){
		foreach my $loc(sort keys %{ $introns{$chr} }){
			my ($ss,$ee) = split(/\t/,$loc);
			my @values;
			my $pass = 0;
			foreach my $pos(sort keys %{ $introns{$chr}{$loc} }){
				my $v = 0;
				if(!$introns{$chr}{$loc}{$pos}){
					$v = 0;
				}else{
					$v = $introns{$chr}{$loc}{$pos};
					$pass++;
				}	
				push(@values,$v);
			}
			next if($pass < 3 || scalar @values == 0);
			my $median = sprintf("%d",median(@values));
			print OUT "$chr\t$loc\t" . $median ."\n";
    	}
    }
    close(OUT);
}
	my $outfn = $bam;
    $outfn=~s/\.Aligned(.*)//;
    $outfn=~s/\.bam(.*)//;
    $outfn .= ".IR.out.tab";
	system("cat \_wososatmp\/" . $bam . ".chr\*.txt > " . $outfn);
	system("rm \_wososatmp\/" . $bam . ".chr\*.txt");
    my $stoptime = time;
    my $mins = ($stoptime-$starttime) / 60;
    print "Spent $mins mins\n";
    
