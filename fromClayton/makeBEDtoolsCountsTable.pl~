$outfile = $ARGV[1] . "bedtools.all.counts.txt";
open(OUT, ">$outfile");
opendir (DIR, $ARGV[0]);
@files = grep /bedtools.counts/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.bedtools\.counts/);
	$filename_hash{$filename} = 1;
	$infile = $ARGV[0] . $file;
	open(IN, $infile);
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
		$geneexon = $line[3];
		$counts = $line[6];
		@geneexon = split(":", $geneexon);
		$gene = $geneexon[0];
		if ($filecount == 1) {
			$exon_start = $line[1];
			$exon_end = $line[2];
			$exon_length = $exon_end - $exon_start + 1;
			$length{$gene} += $exon_length;
		}
		$gene_hash{$gene} = 1;
		$counts{$filename}{$gene} += $counts;
	}
	close(IN);
} 
foreach $i (sort {$a cmp $b} keys %filename_hash) {
	$i =~ s/-/\./g;
	print OUT "\t$i";
}
print OUT "\tLength\n";
foreach $j (sort {$a cmp $b} keys %gene_hash) {
	if ($j eq '') {
		next;
	}
	print OUT "$j";
	foreach $i (sort {$a cmp $b} keys %filename_hash) {
		print OUT "\t$counts{$i}{$j}";
	}
	print OUT "\t$length{$j}\n";
}
close(OUT);
	
	
