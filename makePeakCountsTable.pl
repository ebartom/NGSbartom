use strict;
my $inputDirectory = $ARGV[0];
my $outputDirectory = $ARGV[1];
my $prefix = $ARGV[2];
my $outfile = "$outputDirectory/$prefix.all.counts.txt";
open(OUT, ">$outfile");
opendir (DIR, $inputDirectory);
print STDERR "Looking for $prefix.*.counts.bed in $inputDirectory\n";
my @files = grep /$prefix.*.counts.bed/, readdir DIR;
close DIR;
print STDERR "Found @files\n";
my $filecount = 0;
my (%filename_hash,%gene_hash);
my (%length,%info,%strands,%starts,%ends,%chr,%counts,%hash);
foreach my $file (@files) {
	$filecount++;
#	print STDERR "File: $file\n";
	my ($filename) = ($file =~ /$prefix\.(.+)\.counts\.bed/);
#	print STDERR "F:$filename\n";
#	my $filename = $file;
	$filename_hash{$filename} = 1;
	my $infile = "$inputDirectory/$file";
	open(IN, $infile);
	while (<IN>) {
	    my $line = $_;
	    chomp $line;
#	    print "$line\n";
	    my @line = split("\t", $line);
	    my $chr = $line[0];
	    my $name = $line[3];
	    my $strand = "+";
	    my $counts = $line[4];
	    my $start = $line[1];
	    my $end = $line[2];
	    my $length = $end - $start + 1;	
	    $length{$name} = $length;
	    $hash{$name} = 1;
	    $counts{$filename}{$name} = $counts;
	    $strands{$name} = $strand;
	    $starts{$name} = $start;
	    $ends{$name} = $end;
	    $chr{$name} = $chr;
	}
	close(IN);
	
}
print OUT "Chr\tStart\tEnd\tLength\tStrand";
foreach my $i (sort {$a cmp $b} keys %filename_hash) {
#    $i =~ s/-/\./g;
    print OUT "\t$i";
}
print OUT "\n";


foreach my $j (sort {$starts{$a} <=> $starts{$b} } keys %hash) {
    if ($j eq '') {
	next;
    }
    print OUT "$j";
    print OUT "\t$chr{$j}\t$starts{$j}\t$ends{$j}\t$length{$j}\t$strands{$j}";
    foreach my $i (sort {$a cmp $b} keys %filename_hash) {
#	print STDERR "i,j $i,$j\t$counts{$i}{$j}\n";
	print OUT "\t$counts{$i}{$j}";
    }
    print OUT "\n";
}
close(OUT);

	
