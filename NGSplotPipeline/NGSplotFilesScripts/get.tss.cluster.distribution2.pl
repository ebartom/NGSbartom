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
	if ($strand eq '+') {
		$tss = $start;
	} else {
		$tss = $stop;
	}
	push @{$tss{$chr}}, $tss;
}
close(IN);
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	@{$sorted_tss{$i}} = sort {$a <=> $b} @{$tss{$i}}
}
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	$chr = $i;
	foreach $j (@{$sorted_tss{$i}}) {
		foreach $k (@{$sorted_tss{$i}}) {
			$l = abs($k - $j);
			if ($l == 0 || $l > 100000) {
				next;
			}
			if ($l <= 1000) {
				$count1k++;
			} elsif ($l <= 2000) {
				$count2k++;
			} elsif ($l <= 3000) {
				$count3k++;
			} elsif ($l <= 4000) {
				$count4k++;
			} elsif ($l <= 5000) {
				$count5k++;
			} elsif ($l <= 6000) {
				$count6k++;
			} elsif ($l <= 7000) {
				$count7k++;
			} elsif ($l <= 8000) {
				$count8k++;
			} elsif ($l <= 9000) {
				$count9k++;
			} elsif ($l <= 10000) {
				$count10k++;
			} elsif ($l <= 15000) {
				$count15k++;
			} elsif ($l <= 20000) {
				$count20k++;
			} elsif ($l <= 25000) {
				$count25k++;
			} elsif ($l <= 30000) {
				$count30k++;
			} elsif ($l <= 35000) {
				$count35k++;
			} elsif ($l <= 40000) {
				$count40k++;
			} elsif ($l <= 45000) {
				$count45k++;
			} elsif ($l <= 50000) {
				$count50k++;
			} elsif ($l <= 55000) {
				$count55k++;
			} elsif ($l <= 60000) {
				$count60k++;
			} elsif ($l <= 65000) {
				$count65k++;
			} elsif ($l <= 70000) {
				$count70k++;
			} elsif ($l <= 75000) {
				$count75k++;
			} elsif ($l <= 80000) {
				$count80k++;
			} elsif ($l <= 85000) {
				$count85k++;
			} elsif ($l <= 90000) {
				$count90k++;
			} elsif ($l <= 95000) {
				$count95k++;
			} elsif ($l <= 100000) {
				$count100k++;
			}
		}
	}
}
open(OUT, ">dm3.ensembl.tss.cluster.distribution2.txt");
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
