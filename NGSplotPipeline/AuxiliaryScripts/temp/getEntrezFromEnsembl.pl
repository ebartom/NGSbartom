$entrez{"hg19"} = "/projects/b1025/tools/homer/data/accession/human.description";
$entrez{"hg18"} = "/projects/b1025/tools/homer/data/accession/human.description";
$entrez{"dm3"} = "/projects/b1025/tools/homer/data/accession/fly.description";
$entrez{"mm10"} = "/projects/b1025/tools/homer/data/accession/mouse.description";
$entrez{"mm9"} = "/projects/b1025/tools/homer/data/accession/mouse.description";
$entrez{"sacCer3"} = "/projects/b1025/tools/homer/data/accession/yeast.description";
$assembly = $ARGV[1];
open(IN, $entrez{$assembly});
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$gid = $line[3];
	$symbol{$gid} = $line[4];
	$entrez{$gid} = $line[0];
	print "$line[0]\t$line[3]\t$line[4]\n";
}
close(IN);
opendir(DIR, $ARGV[0]);
@files = grep /txt/, readdir DIR;
close DIR;
foreach $file (@files) {	
	($filename) = ($file =~ /(\S+)\.txt/);
	$entrez = $ARGV[0] . $filename . ".entrez.txt";
	open(OUT, ">$entrez");
	$infile = $ARGV[0] . $file;
	open(IN2, $infile);
	$linecount = 0;
	while ($line = <IN2>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		$gene = $line[0];
		if ($entrez{$gene} eq '') {
			next;
		}
		$count{$filename}{$entrez{$gene}}++;
		if ($count{$filename}{$entrez{$gene}} > 1) {
			next;
		}
		print OUT "$entrez{$gene}\n";
	}
	close(IN2);
	close(OUT);
}
