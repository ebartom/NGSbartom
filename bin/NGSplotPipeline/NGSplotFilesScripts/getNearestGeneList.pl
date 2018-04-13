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
($basename) = ($ARGV[1] =~ /(\S+)\.bed/);
$outfile = $basename . ".nearestGeneList.txt";
open(OUT, ">$outfile");
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
			$closest_gene3 = $tid{$chr}{$i};
    		}
    	}
	$k = $closest_gene;
	$k2 = $closest_gene2;
	$k3 = $closest_gene3;
	$uid = $chr . "_" . $center;
    	#print OUT "$uid\t$k\t$k2\t$k3\n";
    	print OUT "$k3\n";
}
close(IN);
close(OUT);
$clusters = $ARGV[2];
foreach $i (1 .. $clusters) {
	$infile = $basename . ".cluster." . $i . ".bed";
	$outfile = $basename . ".cluster." . $i . ".nearestGeneList.txt";
	open(OUT, ">$outfile");
	open(IN, $infile);
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
		$chr = $line[0];
		$start = $line[1];
		$stop = $line[2];
		$center = int(($stop - $start) / 2) + $start;
    		$max = 300000000;
    		$closest_gene = "";
		$closest_gene2 = "";
		$closest_gene3 = "";
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
		print OUT "$uid\t$k\t$k2\t$k3\n";
	}
	close(IN);
	close(OUT);
}
