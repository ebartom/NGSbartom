opendir(DIR, $ARGV[0]);
@files = grep /hg18.ensembl.genebody/, readdir DIR;
close DIR;
foreach $file (@files) {
	($filename) = ($file =~ /(\S+)\.txt/);
	$outfile = $filename . ".bed";
	open(OUT, ">$outfile");
	$infile = $ARGV[0] . $file;
	open(IN, $infile);
	$linecount = 0; 
	while($line = <IN>) {
		$linecount++;
		if ($linecount == 1) {
			next;
		}
		chomp $line;
		@line = split("\t", $line);
		$chr = $line[1];
		$start = $line[2];
		$stop = $line[3];
		$gid = $line[4];
		$tid = $line[6];
		$strand = $line[7];
		print OUT "$chr\t$start\t$stop\t$tid\t$gid\t$strand\n";
	}
	close(OUT);
	close(IN);
}
