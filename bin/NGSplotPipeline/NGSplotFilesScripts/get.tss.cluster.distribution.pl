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
	$tss_hash{$chr}{$tss} = 1;
}
close(IN);
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	@{$sorted_tss{$i}} = sort {$a <=> $b} @{$tss{$i}}
}
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	$chr = $i;
	foreach $j (@{$sorted_tss{$i}}) {
		foreach $k (-100000 ... 100000) {
			$l = $k + $j;
			if ($tss_hash{$chr}{$l} == 1) {
				$count100k++;
			}
			if (abs{$k} <= 95000 && $tss_hash{$chr}{$l} == 1) {
				$count95k++;
			}
			if (abs{$k} <= 90000 && $tss_hash{$chr}{$l} == 1) {
				$count90k++;
			}
			if (abs{$k} <= 85000 && $tss_hash{$chr}{$l} == 1) {
				$count85k++;
			}
			if (abs{$k} <= 80000 && $tss_hash{$chr}{$l} == 1) {
				$count80k++;
			}
			if (abs{$k} <= 75000 && $tss_hash{$chr}{$l} == 1) {
				$count75k++;
			}
			if (abs{$k} <= 70000 && $tss_hash{$chr}{$l} == 1) {
				$count70k++;
			}
			if (abs{$k} <= 65000 && $tss_hash{$chr}{$l} == 1) {
				$count65k++;
			}
			if (abs{$k} <= 60000 && $tss_hash{$chr}{$l} == 1) {
				$count60k++;
			}
			if (abs{$k} <= 55000 && $tss_hash{$chr}{$l} == 1) {
				$count55k++;
			}
			if (abs{$k} <= 50000 && $tss_hash{$chr}{$l} == 1) {
				$count50k++;
			}
			if (abs{$k} <= 45000 && $tss_hash{$chr}{$l} == 1) {
				$count45k++;
			}
			if (abs{$k} <= 40000 && $tss_hash{$chr}{$l} == 1) {
				$count40k++;
			}
			if (abs{$k} <= 35000 && $tss_hash{$chr}{$l} == 1) {
				$count35k++;
			}
			if (abs{$k} <= 30000 && $tss_hash{$chr}{$l} == 1) {
				$count30k++;
			}
			if (abs{$k} <= 25000 && $tss_hash{$chr}{$l} == 1) {
				$count25k++;
			}
			if (abs{$k} <= 20000 && $tss_hash{$chr}{$l} == 1) {
				$count20k++;
			}
			if (abs{$k} <= 15000 && $tss_hash{$chr}{$l} == 1) {
				$count15k++;
			}
			if (abs{$k} <= 10000 && $tss_hash{$chr}{$l} == 1) {
				$count10k++;
			}
			if (abs{$k} <= 9000 && $tss_hash{$chr}{$l} == 1) {
				$count9k++;
			}
			if (abs{$k} <= 8000 && $tss_hash{$chr}{$l} == 1) {
				$count8k++;
			}
			if (abs{$k} <= 7000 && $tss_hash{$chr}{$l} == 1) {
				$count7k++;
			}
			if (abs{$k} <= 6000 && $tss_hash{$chr}{$l} == 1) {
				$count6k++;
			}
			if (abs{$k} <= 5000 && $tss_hash{$chr}{$l} == 1) {
				$count5k++;
			}
			if (abs{$k} <= 4000 && $tss_hash{$chr}{$l} == 1) {
				$count4k++;
			}
			if (abs{$k} <= 3000 && $tss_hash{$chr}{$l} == 1) {
				$count3k++;
			}
			if (abs{$k} <= 2000 && $tss_hash{$chr}{$l} == 1) {
				$count2k++;
			}
			if (abs{$k} <= 1000 && $tss_hash{$chr}{$l} == 1) {
				$count1k++;
			}
		}
	}
}
open(OUT, ">dm3.ensembl.tss.cluster.distribution.txt");
print OUT "1000\t$count1k\n";
print OUT "2000\t$count2k\n";
print OUT "3000\t$count3k\n";
print OUT "4000\t$count4k\n";
print OUT "5000\t$count5k\n";
print OUT "6000\t$count6k\n";
print OUT "7000\t$count7k\n";
print OUT "8000\t$count8k\n";
print OUT "9000\t$count9k\n";
print OUT "10000\t$count10k\n";
print OUT "15000\t$count15k\n";
print OUT "20000\t$count20k\n";
print OUT "25000\t$count25k\n";
print OUT "30000\t$count30k\n";
print OUT "35000\t$count35k\n";
print OUT "40000\t$count40k\n";
print OUT "45000\t$count45k\n";
print OUT "50000\t$count50k\n";
print OUT "55000\t$count55k\n";
print OUT "60000\t$count60k\n";
print OUT "65000\t$count65k\n";
print OUT "70000\t$count70k\n";
print OUT "75000\t$count75k\n";
print OUT "80000\t$count80k\n";
print OUT "85000\t$count85k\n";
print OUT "90000\t$count90k\n";
print OUT "95000\t$count95k\n";
print OUT "100000\t$count100k\n";
close(OUT);
