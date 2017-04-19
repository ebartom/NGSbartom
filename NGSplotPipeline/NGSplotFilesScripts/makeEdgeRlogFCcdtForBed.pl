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
$edgeRpv = $ARGV[4];
open(IN, $ARGV[1]);
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	$infile = $line[1];
	open(IN2, $infile);
	$linecount = 0;
	while ($line = <IN2>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				}
				if ($line[$i] =~ /adj/) {
					$adjp = $i;
				}
			}
			next;
		}
		$gene = $line[0];
		$logfc{$name}{$gene} = $line[$logfc];
		if ($logfc{$name}{$gene} < 0 && $line[$adjp] < $edgeRpv) {
			$down{$name}{$gene} = 1;
		}
		if ($logfc{$name}{$gene} > 0 && $line[$adjp] < $edgeRpv) {
			$up{$name}{$gene} = 1;
		}
	}
	close(IN2);
}
close(IN);
open(IN, $ARGV[2]);
($basename) = ($ARGV[2] =~ /(\S+)\.bed/);
$outfile = $basename . ".edgeR.logFC.cdt";
$outfile2 = $basename . ".edgeR.logFC.geneNames.cdt";
open(OUT, ">$outfile");
print OUT "UID\tNAME";
foreach $i (@names) {
	print OUT "\t$i";
}
print OUT "\n";
open(OUT2, ">$outfile2");
print OUT2 "UID\tNAME";
foreach $i (@names) {
	print OUT2 "\t$i";
}
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
	foreach $j (@names) {
		if ($logfc{$j}{$k} eq '') {
			$logfc{$j}{$k} = 0;
		}
		#use k, never k2 for geneID
		print OUT "\t$logfc{$j}{$k}";
		print OUT2 "\t$logfc{$j}{$k}";
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
	$outfile = $basename . ".cluster." . $i . ".edgeR.logFC.cdt";
	$outfile2 = $basename . ".cluster." . $i . ".edgeR.logFC.geneNames.cdt";
	foreach $j (@names) {
		$outfile3 = $basename . ".cluster." . $i . "." . $j . ".edgeR.down." . $edgeRpv . ".bed";
		$filehandledown{$i}{$j} = "down." . $i . "." . $j;
		open($filehandledown{$i}{$j}, ">$outfile3");
	}
	foreach $j (@names) {
		$outfile4 = $basename . ".cluster." . $i . "." . $j . ".edgeR.up." . $edgeRpv . ".bed";
		$filehandleup{$i}{$j} = "up." . $i . "." . $j;
		open($filehandleup{$i}{$j}, ">$outfile4");
	}
	open(IN, $infile);
	open(OUT, ">$outfile");
	print OUT "UID\tNAME";
	foreach $e (@names) {
		print OUT "\t$e";
	}
	print OUT "\n";
	open(OUT2, ">$outfile2");
	print OUT2 "UID\tNAME";
	foreach $e (@names) {
		print OUT2 "\t$e";
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
			if ($logfc{$d}{$k} eq '') {
				$logfc{$d}{$k} = 0;
			}
			#use k, never k2 for geneID
			print OUT "\t$logfc{$d}{$k}";
			print OUT2 "\t$logfc{$d}{$k}";
			if ($down{$d}{$k} == 1) {
				$name = $chr . "_" . $center;
				print {$filehandledown{$i}{$d}} "$chr\t$start\t$stop\t$name\t$k\t$strand{$k3}\n";
			}
			if ($up{$d}{$k} == 1) {
				$name = $chr . "_" . $center;
				print {$filehandleup{$i}{$d}} "$chr\t$start\t$stop\t$name\t$k\t$strand{$k3}\n";
			}	
		}
		print OUT "\n";
		print OUT2 "\n";
	}
	close(IN);
	close(OUT);
	close(OUT2);
}
