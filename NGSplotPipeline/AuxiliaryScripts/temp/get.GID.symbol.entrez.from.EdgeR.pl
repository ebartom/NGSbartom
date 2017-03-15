$entrez{"hg19"} = "/projects/b1025/tools/homer/data/accession/human.description";
$entrez{"hg18"} = "/projects/b1025/tools/homer/data/accession/human.description";
$entrez{"dm3"} = "/projects/b1025/tools/homer/data/accession/fly.description";
$entrez{"mm10"} = "/projects/b1025/tools/homer/data/accession/mouse.description";
$entrez{"mm9"} = "/projects/b1025/tools/homer/data/accession/mouse.description";
$entrez{"sacCer3"} = "/projects/b1025/tools/homer/data/accession/yeast.description";
$assembly = $ARGV[2];
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
open(IN, $ARGV[0]);
$edgeRpv = $ARGV[1];
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	($file) = ($ARGV[0] =~ /(\S+)\.txt/);
	$outfiledown = $file . "." . $name . ".edgeR.down." . $edgeRpv . ".all.txt";
	$outfiledownentrez = $file . "." . $name . ".edgeR.down." . $edgeRpv . ".all.entrez.txt";
	$outfileup = $file . "." . $name . ".edgeR.up." . $edgeRpv . ".all.txt";
	$outfileupentrez = $file . "." . $name . ".edgeR.up." . $edgeRpv . ".all.entrez.txt";  
	open(DOWN, ">$outfiledown");
	open(UP, ">$outfileup");
	open(DOWNE, ">$outfiledownentrez");
	open(UPE, ">$outfileupentrez");
	$infile = $line[1];
	open(IN2, $infile);
	$linecount = 0;
	while ($line = <IN2>) {
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
				if ($line[$i] =~ /gene/) {
					$symbol = $i;
				}
			}
			next;
		}
		$gene = $line[0];
		$logfc{$name}{$gene} = $line[$logfc];
		$genename = $line[$symbol];
		if ($logfc{$name}{$gene} < 0 && $line[$adjp] < $edgeRpv) {
			$down{$name}{$gene} = 1;
			print DOWN "$gene	$genename	$entrez{$gene}\n";
			if ($entrez{$gene} eq '') {
			} else {
			print DOWNE "$entrez{$gene}\n";
			}
		}
		if ($logfc{$name}{$gene} > 0 && $line[$adjp] < $edgeRpv) {
			$up{$name}{$gene} = 1;
			print UP "$gene	$genename	$entrez{$gene}\n";
			if ($entrez{$gene} eq '') {
			} else {
			print UPE "$entrez{$gene}\n";
			}
		}
	}
	close(IN2);
}
close(IN);
close(DOWN);
close(DOWNE);
close(UP);
close(UPE);
