open(IN, $ARGV[0]);
$linecount = 0;
while ($line = <IN>) {
	$linecount++;
	if ($linecount == 1) {
		next;
	}
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$chr = $line[1];
	$chr_hash{$chr} = 1;
	$start = $line[2];
	$stop = $line[3];
	$tid = $line[6];
	$gid = $line[4];
	$strand = $line[7];
	$gene = $line[5];
	if ($strand eq '+') {
		$tss = $start;
	} else {
		$tss = $stop;
	}
	push @{$tss{$chr}}, $tss;
	$tid{$chr}{$tss} = $tid; 
	$gid{$tid} = $gid;
	$gene{$tid} = $gene;
	$strand{$tid} = $strand;
}
close(IN);
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	@{$sorted_tss{$i}} = sort {$a <=> $b} @{$tss{$i}}
}
open(IN, $ARGV[1]);
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	$infile = $line[1];
	print "$name	$infile\n";
	open(IN2, $infile);
	$linecount = 0;
	$logfc = 0;
	while ($line = <IN2>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				} elsif ($logfc == 0 && $i > 5) {
					push @{$samples{$name}}, $line[$i];
				} else {
				}
				#print "$logfc";
				#foreach $j (@{$samples{$name}}) {
				#	print "\t$j";
				#}
				#print "\n";
			}
			next;
		}
		$gene = $line[0];
		foreach $i (6 .. $logfc - 1) {
			push @{$logCPM{$name}{$gene}}, $line[$i];
		}
	}
	close(IN2);
}

close(IN);
open(IN, $ARGV[2]);
($basename) = ($ARGV[2] =~ /(\S+)\.bed/);
$outfile = $basename . ".edgeR.logCPM.cdt";
$outfile2 = $basename . ".edgeR.logCPMgeneNames.cdt";
open(OUT, ">$outfile");
print OUT "UID\tNAME";
open(OUT2, ">$outfile2");
print OUT2 "UID\tNAME";
foreach $i (@names) {
	foreach $j (@{$samples{$i}}) {
		$temp = $i . "." . $j;
		print "$temp\n";
        	print OUT "\t$temp";
		print OUT2 "\t$temp";
	}
}
print OUT "\n";
print OUT2 "\n";
while($line = <IN>){;
	chomp $line;
	@line = split("\t", $line);
	$chr = $line[0];
	$start = $line[1];
	$stop = $line[2];
	$center = int(($stop - $start) / 2) + $start;
   	$max = 300000000;
    $closest_gene = "";
    $closest_gene2 = "";
    foreach $i (@{$sorted_tss{$chr}}) {
    	$distance = abs($center - $i);
    	if ($distance < $max) {
    		$max = $distance;
    		$closest_gene = $gid{$tid{$chr}{$i}};
    		$closest_gene2 = $gene{$tid{$chr}{$i}};
    	}
    }
	$k = $closest_gene;
	$k2 = $closest_gene2;
	$uid = $chr . "_" . $center;
    print OUT "$uid\t$k";
    print OUT2 "$uid\t$k2";
	foreach $d (@names) {
			#use k, never k2 for geneID
			@array = @{$samples{$d}};
			foreach $y (0 .. $#array) {
				$x = $logCPM{$d}{$k}[$y];
				if ($x eq '') {
					$x = 0;
				}
				print OUT "\t$x";
				print OUT2 "\t$x";
			}
	}
	print OUT "\n";
	print OUT2 "\n";
}
close(IN);
close(OUT);
close(OUT2);
$clusters = $ARGV[3];
if ($clusters == 0 || $clusters eq '') {
	exit;
}
foreach $i (1 .. $clusters) {
	$infile = $basename . ".cluster." . $i . ".bed";
	$outfile = $basename . ".cluster." . $i . ".edgeR.logCPM.cdt";
	$outfile2 = $basename . ".cluster." . $i . ".edgeR.logCPM.geneNames.cdt";
	open(IN, $infile);
	open(OUT, ">$outfile");
	print OUT "UID\tNAME";
	foreach $e (@names) {
		foreach $s (@{$samples{$d}}) {
    		$temp = $e . "." . $s;
			print OUT "\t$temp";
		}
	}
	print OUT "\n";
	open(OUT2, ">$outfile2");
	print OUT2 "UID\tNAME";
	foreach $e (@names) {
		foreach $s (@{$samples{$d}}) {
    		$temp = $e . "." . $s;
			print OUT2 "\t$temp";
		}
	}
	print OUT2 "\n";
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
		$chr = $line[0];
		$start = $line[1];
		$stop = $line[2];
		$center = int(($stop - $start) / 2) + $start;
    		$max = 300000000;
    		$closest_gene = "";
    		foreach $c (@{$sorted_tss{$chr}}) {
    			$distance = abs($center - $c);
    			if ($distance < $max) {
    				$max = $distance;
    				$closest_gene = $gid{$tid{$chr}{$c}};
    				$closest_gene2 = $gene{$tid{$chr}{$c}};
    				$closest_gene3 = $tid{$chr}{$c};
    			}
   		}
		$k = $closest_gene;
		$k2 = $closest_gene2;
		$k3 = $closest_gene3;
		$uid = $chr . "_" . $center;
		print OUT "$uid\t$k";
		print OUT2 "$uid\t$k2";
		foreach $d (@names) {
			#use k, never k2 for geneID
			@array = @{$samples{$d}};
			foreach $y (0 .. $#array) {
				$x = $logCPM{$d}{$k}[$y];
				if ($x == '') {
					$x = 0;
				}
				print OUT "\t$x";
				print OUT2 "\t$x";
			}
		}
		print OUT "\n";
		print OUT2 "\n";
	}
	close(IN);
	close(OUT);
	close(OUT2);
}
