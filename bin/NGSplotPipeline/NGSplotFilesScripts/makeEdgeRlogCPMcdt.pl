open(IN, $ARGV[0]);
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
open(IN, $ARGV[1]);
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	$infile = $line[1];
	#print "$name\t$infile\n";
	open(IN2, $infile);
	$linecount = 0;
	$logfc eq ''
	while ($line = <IN2>) {
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				} elsif ($logfc eq '' & $i > 5) {
					push @{$samples{$name}}, $line[$i];
				}
			}
			next;
		}
		$gene = $line[0];
		foreach $i (6 .. $logfc - 1) {
			push @{$logCPM{$name}{$gene}}, $line[$i];
		}
	}
	close(IN2);
}
close(IN);
($basename) = ($ARGV[2] =~ /(\S+)\.txt/);
$outfile = $basename . ".edgeR.logCPM.cdt";
$outfile2 = $basename . ".edgeR.logCPM.geneNames.cdt";
open(OUT, ">$outfile");
print OUT "UID\tNAME";
open(OUT2, ">$outfile2");
print OUT2 "UID\tNAME";
foreach $i (@names) {
	foreach $j (@{$samples{$name}}) {
		$temp = $i . "." . $j;
        print OUT "\t$temp";
		print OUT2 "\t$temp";
	}
}
print OUT1 "\n";
print OUT2 "\n";
open(IN, $ARGV[2]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[0];
	$i = $tid;
	$j = $gene{$i};
	print OUT "$i\t$gid{$i}";
	print OUT2 "$i\t$j";
	foreach $j (@names) {
		foreach $x (@{$logCPM{$j}{$gid{$i}}}) {
			if ($x eq '') {
				$x = 0;
			}
			print OUT "\t$x";
			print OUT2 "\t$x";
		}
	}
	print OUT "\n";
	print OUT2 "\n";
}
close(IN);
close(OUT);
$clusters = $ARGV[3];
foreach $i (1 .. $clusters) {
	$infile = $basename . ".cluster." . $i . ".txt";
	$outfile = $basename . ".cluster." . $i . ".edgeR.logCPM.cdt";
	$outfile2 = $basename . ".cluster." . $i . ".edgeR.logCPM.geneNames.cdt";
	open(IN, $infile);
	open(OUT, ">$outfile");
	print OUT "UID\tNAME";
	foreach $i (@names) {
		foreach $s (@{$samples{$i}}) {
    			$temp = $i . "." . $s;
			print OUT "\t$temp";
		}
	}
	print OUT "\n";
	open(OUT2, ">$outfile2");
	print OUT2 "UID\tNAME";
	foreach $i (@names) {
		foreach $s (@{$samples{$i}}) {
            $temp = $i . "." . $s;
            print OUT2 "\t$temp";
        }
	}
	print OUT2 "\n";
	while ($line = <IN>) {
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		$tid = $line[0];
		$i = $tid;
		$j = $gene{$i};
		print OUT "$i\t$gid{$i}";
		print OUT2 "$i\t$j";
		foreach $d (@names) {
			foreach $x (@{$logCPM{$d}{$gid{$i}}}) {
				if ($x eq '') {
					$x = 0;
				}
				print OUT "\t$x";
				print OUT2 "\t$x";
			}		
		}
		print OUT "\n";
		print OUT2 "\n";
	}
	close(IN);
	close(OUT);
	close(OUT2);
}
