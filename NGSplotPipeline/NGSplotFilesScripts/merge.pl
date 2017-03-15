open(OUT, ">merge.sh");
print OUT "module load samtools\n";
open(IN, $ARGV[0]);
while($line=<IN>) {
	chomp $line;
	@line = split("\t", $line);
	$outfile = $line[0] . ".bam";
	$file_hash{$line[0]} = 1;
	$infile1 = $line[1];
	$infile2 = $line[2];
	print OUT "samtools merge $outfile $infile1 $infile2 \&\n";
}
print OUT "wait\n";
foreach $i (sort {$a cmp $b} keys %file_hash) {
	$infile = $i . ".bam";
	$outfile = $i . ".sorted";
	print OUT "samtools sort $infile $outfile \&\n";
}
print OUT "wait\n";
foreach $i (sort {$a cmp $b} keys %file_hash) {
	$infile = $i . ".sorted.bam";
	print OUT "samtools index $infile \&\n";
}
print OUT "wait\n";


	