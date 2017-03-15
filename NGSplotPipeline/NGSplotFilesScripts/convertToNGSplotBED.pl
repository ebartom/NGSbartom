open(OUT, ">$ARGV[1]");
open(IN, $ARGV[0]);
while($line = <IN>){;
	chomp $line;
	@line = split("\t", $line);
	$chr = $line[0];
	$start = $line[1];
	$stop = $line[2];
	$width = $stop - $start + 1;
	$center = int(($stop - $start) / 2) + $start;
    	#$center1 = $center + 1;
	$name = $chr . "_" . $center;
	$strand = $line[5];
	if ($strand eq '') {
		$strand = '+';
	}
	print OUT "$chr\t$start\t$stop\t$name\t$width\t$strand\n";
}
print "\t$linecount\n";
close(IN);
close(OUT);
