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
	$width{$name} = $width;
	$chr{$name} = $chr;
	$start{$name} = $start;
	$stop{$name} = $stop;
	#$center{$name} = $center;
	#$center1{$name} = $center1;
	$strand{$name} = $strand;
}
foreach $i (sort {$width{$b} <=> $width{$a}} keys %width) {
	#print OUT "$chr{$i}\t$center{$i}\t$center1{$i}\t$i\t$width{$i}\t$strand{$i}\n";
	print OUT "$chr{$i}\t$start{$i}\t$stop{$i}\t$i\t$width{$i}\t$strand{$i}\n";
}
print "\t$linecount\n";
close(IN);
close(OUT);
