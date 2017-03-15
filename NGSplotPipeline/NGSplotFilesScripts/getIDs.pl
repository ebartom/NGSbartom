open(OUT, ">$ARGV[1]");
open(IN, $ARGV[0]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	@go = split(":", $line[0]);
	$id = $go[0];
	print OUT "$id\n";
}
close(IN);
close(OUT);
