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
        $id = $line[3];
        $bed{$id} = $line;;
}
close(IN);
open(IN, $ARGV[2]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$id = $line[0];
	push @ids, $id;
}
close(IN);
$count = 0;
($file) = ($ARGV[2] =~ /(\S+)\.txt/);
$outfile = $file . ".bed";
open(ALL, ">$outfile");
foreach $i (sort {$a <=> $b} keys %cluster_count) {
	($file) = ($ARGV[2] =~ /(\S+)\.txt/);
	$outfile = $file . ".cluster." . $i . ".bed"; 
	open(OUT, ">$outfile");
	foreach $j (1 .. $cluster_count{$i}) {
		print OUT "$bed{$ids[$count]}\n";
		print ALL "$bed{$ids[$count]}\n";
		$count++;
	}
	close(OUT);
}		
