$outfile = $ARGV[2] . "htseq.all.counts.txt";
open(OUT, ">$outfile");
open(IN, $ARGV[1]);
while ($line = <IN>) {
	chomp $line;
	@line = split("\t", $line);
	$feature = $line[2];
	if ($feature eq 'exon') {
		$exon_start = $line[3];
		$exon_end = $line[4];
		$exon_length = $exon_end - $exon_start + 1;
		$info = $line[8];
		@info = split(";", $info);
		$id = $info[0];
		$id =~ s/"//g;
		@id = split(" ", $id);
		$gene = $id[1];
		#print "$gene	$exon_length\n";
		$length{$gene} += $exon_length;
	}
}
close(IN);
opendir (DIR, $ARGV[0]);
@files = grep /htseq.counts/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.htseq\.counts/);
	$filename_hash{$filename} = 1;
	$infile = $ARGV[0] . $file;
	open(IN, $infile);
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
		$gene = $line[0];
		$counts = $line[1];
		if ($gene eq '__no_feature' || $gene eq '__ambiguous' || $gene eq '__too_low_aQual' || $gene eq '__not_aligned' || $gene eq '__alignment_not_unique') {
			next;
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
	
