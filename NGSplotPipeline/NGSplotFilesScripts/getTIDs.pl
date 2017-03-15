open(IN, $ARGV[2]);
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
}
close(IN);
open(IN, $ARGV[3]);
$linecount = 0;
while ($line = <IN>) {
	$linecount++;
	if ($linecount == 1) {
		next;
	}
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$entrez{$line[3]} = $line[0];
}
close(IN);
open(OUT, ">$ARGV[1]");
open(IN, $ARGV[0]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	@go = split(":", $line[0]);
	$tid = $go[1];
	print OUT "$tid	$gid{$tid}	$gene{$tid}	$entrez{$gid{$tid}}\n";
}
close(IN);
close(OUT);