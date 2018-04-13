use strict;
my $inputDirectory = $ARGV[0];
my $gtfFile = $ARGV[1];
my $outputDirectory = $ARGV[2];

my $outfile = "$outputDirectory/htseq.all.counts.txt";
print STDERR "$outfile\n";
open(OUT, ">$outfile");
open(IN, $gtfFile);
my (%length,%info,%strands,%starts,%ends,%chr);
while (my $line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	my $feature = $line[2];
	if ($feature eq 'exon') {
		my $chr = $line[0];
		my $exon_start = $line[3];
		my $exon_end = $line[4];
		my $exon_length = $exon_end - $exon_start + 1;
		my $strand = $line[6];
		my $info = $line[8];
		my $gene;
#		my @info = split(";", $info);
#		my $id = $info[0];
#		$id =~ s/"//g;
#		my @id = split(" ", $id);
		#		my $gene = $id[1];
		if ($info =~ /gene_id \"([\w\d\-\.\_\)\(]+)\"/){
		    $gene = $1;
		}
#		print STDERR "$info\n$gene\n";
		#print "$gene	$exon_length\n";
		$length{$gene} += $exon_length;
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
}
close(IN);
opendir (DIR, $inputDirectory);
my @files = grep /htseq.counts/, readdir DIR;
close DIR;
my $filecount = 0;
my (%filename_hash,%gene_hash,%counts);
foreach my $file (@files) {
	$filecount++;
	print "$file\n";
	my ($filename) = ($file =~ /(\S+)\.htseq\.counts/);
	$filename_hash{$filename} = 1;
	my $infile = "$inputDirectory/$file";
	open(IN, $infile);
	while (my $line = <IN>) {
		chomp $line;
		my @line = split("\t", $line);
		my $gene = $line[0];
		my $counts = $line[1];
		if ($gene eq '__no_feature' || $gene eq '__ambiguous' || $gene eq '__too_low_aQual' || $gene eq '__not_aligned' || $gene eq '__alignment_not_unique') {
			next;
		}
		$gene_hash{$gene} = 1;
		$counts{$filename}{$gene} += $counts;
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
	
