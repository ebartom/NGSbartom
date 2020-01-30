$edgeRpv = $ARGV[1];
opendir(DIR, $ARGV[0]);
@files = grep /edgeR.txt/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.edgeR\.txt/);
	$outfiledown = $ARGV[0] . $file . "." . $name . ".edgeR.down." . $edgeRpv . ".all.txt";
	$outfileup = $ARGV[0] . $file . "." . $name . ".edgeR.up." . $edgeRpv . ".all.txt";
	open(DOWN, ">$outfiledown");
	open(UP, ">$outfileup");
	$infile = $ARGV[0] . $file;
	open(IN, $infile);
	$linecount = 0;
	while ($line = <IN>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				}
				if ($line[$i] =~ /adj/) {
					$adjp = $i;
				}
			}
			next;
		}
		$gene = $line[0];
		$logfc{$name}{$gene} = $line[$logfc];
		if ($logfc{$name}{$gene} < 0 && $line[$adjp] < $edgeRpv) {
			$down{$name}{$gene} = 1;
			print DOWN "$gene\n";
		}
		if ($logfc{$name}{$gene} > 0 && $line[$adjp] < $edgeRpv) {
			$up{$name}{$gene} = 1;
			print UP "$gene\n";
		}
	}
	close(IN);
}