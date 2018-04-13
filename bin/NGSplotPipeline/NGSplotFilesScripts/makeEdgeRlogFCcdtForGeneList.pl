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
	$gid{$line[6]} = $line[4];
	$gene{$line[6]} = $line[5];
	#print "$line[4]	$line[6]	$gid{$line[6]}\n";
}
close(IN);
open(IN, $ARGV[1]);
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	$infile = $line[1];
	open(IN2, $infile);
	$linecount = 0;
	while ($line = <IN2>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				}
			}
			next;
		}
		$gene = $line[0];
		$logfc{$name}{$gene} = $line[$logfc];
		#print "$gene	$line[$logfc]\n";
	}
	close(IN2);
	#print "$name	$linecount\n";
}
close(IN);
($basename) = ($ARGV[2] =~ /(\S+)\.txt/);
$outfile = $basename . ".edgeR.logFC.cdt";
$outfile2 = $basename . ".edgeR.logFC.geneNames.cdt";
open(OUT, ">$outfile");
print OUT "UID\tNAME";
foreach $i (@names) {
	print OUT "\t$i";
}
print OUT "\n";
open(OUT2, ">$outfile2");
print OUT2 "UID\tNAME";
foreach $i (@names) {
	print OUT2 "\t$i";
}
print OUT2 "\n";
open(IN, $ARGV[2]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	#$line =~ s/\r//g;
	@line = split("\t", $line);
	$tid = $line[0];
	$i = $tid;
	$j = $gene{$i};
	print OUT "$i\t$gid{$i}";
	print OUT2 "$i\t$j";
	foreach $j (@names) {
		if ($logfc{$j}{$gid{$i}} eq '') {
			$logfc{$j}{$gid{$i}} = 0;
		}
		print OUT "\t$logfc{$j}{$gid{$i}}";
		print OUT2 "\t$logfc{$j}{$gid{$i}}";
	}
	print OUT "\n";
	print OUT2 "\n";
}
close(IN);
close(OUT);
close(OUT2);
