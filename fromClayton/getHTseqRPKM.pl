$outfile = $ARGV[1] . "htseq.all.rpkm.txt";
open(OUT, ">$outfile");
open(IN, $ARGV[0]);
$linecount = 0;
while ($line = <IN>) {
	$linecount++;
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	if ($linecount == 1) {
		foreach $i (1 .. ($#line - 1)) {
			print OUT "\t$line[$i]";
			$sample{$i} = $line[$i];
		}
		print OUT "\tLength\n";
		next:
	}
	$gene = $line[0];
	$length{$gene} = $line[$#line];
	$gene_hash{$gene} = 1;
	foreach $i (1 .. ($#line - 1)) {
		$counts{$sample{$i}}{$gene} = $line[$i];
		$total{$sample{$i}} += $line[$i];
	}
}
close(IN);
open(IN, $ARGV[0]);
$linecount = 0;
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	if ($linecount == 1) {
		next;
	}
	$gene = $line[0];
	if ($gene eq '') {
                next;
        }
	print OUT "$gene";
	foreach $i (sort {$a <=> $b} keys %sample) {
		$rpkm = $counts{$sample{$i}}{$gene} / ($total{$sample{$i}} / 1000000) / ($length{$gene} / 1000);
		print OUT "\t$rpkm";
	}
	print OUT "\t$length{$gene}\n";
}
close(IN);
close(OUT);
	
