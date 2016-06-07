use strict;
my $inputDirectory = $ARGV[0];
my $outputDirectory = $ARGV[1];
my $outfile = "$outputDirectory/bedtools.all.counts.txt";
open(OUT, ">$outfile");
opendir (DIR, $inputDirectory);
my @files = grep /bedtools.counts/, readdir DIR;
close DIR;
my $filecount = 0;
my (%filename_hash,%gene_hash);
my (%length,%info,%strands,%starts,%ends,%chr,%gene_hash,%counts);
foreach my $file (@files) {
	$filecount++;
	print "$file\n";
	my ($filename) = ($file =~ /(\S+)\.bedtools\.counts/);
	$filename_hash{$filename} = 1;
	my $infile = "$inputDirectory/$file";
	open(IN, $infile);
	while (my $line = <IN>) {
		chomp $line;
		my @line = split("\t", $line);
		my $chr = $line[0];
		my $geneexon = $line[3];
		my $strand = $line[5];
		my $counts = $line[6];
		my @geneexon = split(":", $geneexon);
		my $gene = $geneexon[0];
		my $exon_start = $line[1];
		my $exon_end = $line[2];
		my $exon_length = $exon_end - $exon_start + 1;		
		if ($filecount == 1) {
			$length{$gene} += $exon_length;
		}
		$gene_hash{$gene} = 1;
		$counts{$filename}{$gene} += $counts;
		$strands{$gene} = $strand;
		if (exists($starts{$gene})){
		    if ($starts{$gene} > $exon_start){ 
			$starts{$gene} = $exon_start;
		    }
		} else { $starts{$gene} = $exon_start;}
		$chr{$gene} = $chr;
		if (exists($ends{$gene})){
		    if ($ends{$gene} < $exon_end){ 
			$ends{$gene} = $exon_end;
		    }
		} else { $ends{$gene} = $exon_end;}
		$chr{$gene} = $chr;
		
	}
	close(IN);
} 
print OUT "Chr\tStart\tEnd\tLength\tStrand";
foreach my $i (sort {$a cmp $b} keys %filename_hash) {
	$i =~ s/-/\./g;
	print OUT "\t$i";
}
print OUT "\n";
foreach my $j (sort {$a cmp $b} keys %gene_hash) {
	if ($j eq '') {
		next;
	}
	print OUT "$j";
	print OUT "\t$chr{$j}\t$starts{$j}\t$ends{$j}\t$length{$j}\t$strands{$j}";
	foreach my $i (sort {$a cmp $b} keys %filename_hash) {
		print OUT "\t$counts{$i}{$j}";
	}
	print OUT "\n";
}
close(OUT);
	
	
