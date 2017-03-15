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
			unless ($l == 0) {
				$count{$l}++;
			}
		}
	}
}
open(OUT, ">dm3.ensembl.tss.cluster.distribution4.txt");
for ($i=1000;$i<=1000000;$i+=1000) {
        $sum = 0;
	foreach $j (0 .. 999) {
		$k = $i - $j;
		$sum += $count{$k};
	}
	print OUT "$i\t$sum\n";
}
