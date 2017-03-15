open(IN, $ARGV[0]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$cluster = $line[1];
	$cluster_count{$cluster}++;
}
close(IN);
open(IN, $ARGV[1]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[0];
	push @tids, $tid;
	$gid{$tid} = $line[1];
	$gene{$tid} = $line[2];
	$entrez{$tid} = $line[3];
}
close(IN);
$count = 0;
foreach $i (sort {$a <=> $b} keys %cluster_count) {
	($file) = ($ARGV[1] =~ /(\S+)\.txt/);
	$outfile = $file . ".cluster." . $i . ".txt"; 
	open(OUT, ">$outfile");
	$outfile2 = $file . ".cluster." . $i . ".entrez.txt"; 
	open(OUT2, ">$outfile2");
	foreach $j (1 .. $cluster_count{$i}) {
		print OUT "$tids[$count]	$gid{$tids[$count]}	$gene{$tids[$count]}	$entrez{$tids[$count]}\n";
		unless ($entrez{$tids[$count]} eq '') {
			print OUT2 "$entrez{$tids[$count]}\n";
		}
		$count++;
	}
	close(OUT);
}

		
