use strict;
my $inputDirectory = $ARGV[0];
my $gtfFile = $ARGV[1];
my $outputDirectory = $ARGV[2];
my $mode = $ARGV[3];

my $outfile = "$outputDirectory/rsem.all.counts.txt";
print STDERR "$outfile\n";
open(OUT, ">$outfile");
my ($outfile2);
if ($mode eq "isoforms"){
    $outfile2 = "$outputDirectory/rsem.all.isopct.txt";
}
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
		my @info = split(";", $info);
		#print "Info:@info\n";
		my $id = "";
		if ($mode eq "genes"){
		    $id = $info[0];
		} elsif( $mode eq "isoforms"){
		    $id = $info[1];
		}
		$id =~ s/"//g;
		my @id = split(" ", $id);
		#print "id:@id\n";
		my $gene = $id[1];
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
my @files;
if ($mode eq "genes"){
    @files = grep /genes.results/, readdir DIR;
} elsif ($mode eq "isoforms"){
    @files = grep /isoforms.results/, readdir DIR;
}
close DIR;
my $filecount = 0;
my (%filename_hash,%gene_hash,%counts,%isopct);
foreach my $file (@files) {
	$filecount++;
	print "$file\n";
	my ($filename) = ($file =~ /(\S+)\.$mode\.results/);
	$filename_hash{$filename} = 1;
	my $infile = "$inputDirectory/$file";
	open(IN, $infile);
	while (my $line = <IN>) {
		chomp $line;
		my @line = split("\t", $line);
		my $gene = $line[0];
		my $counts = sprintf("%.0f",$line[5]);
		$gene_hash{$gene} = 1;
		$counts{$filename}{$gene} = $counts;
		if ($mode eq "isoforms"){
		    $isopct{$filename}{$gene} = $line[7];
		}
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

if ($mode eq "isoforms"){
    open(OUT,">$outfile2");
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
	    print OUT "\t$isopct{$i}{$j}";
	}
	print OUT "\n";
    }
    close(OUT);
}
