open(IN, $ARGV[0]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$cluster = $line[1];
	$cluster_count{$cluster}++;
}
close(IN);
open(IN, $ARGV[1]);
while ($line = <IN>) {
	$line =~ s/"//g;
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[0];
	push @tids, $tid;
	$gid{$tid} = $line[1];
	$gene{$tid} = $line[2];
	$entrez{$tid} = $line[3];
	$tid{$line[1]} = $tid;
}
open(IN, $ARGV[2]);
$edgeRpv = $ARGV[3];
while($line = <IN>) {
	chomp $line;
	my @line = split("\t", $line);
	$name = $line[0];
	push @names, $name;
	($file) = ($ARGV[1] =~ /(\S+)\.txt/);
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
			}
			next;
		}
		$gene = $line[0];
		$logfc{$name}{$gene} = $line[$logfc];
		if ($logfc{$name}{$gene} < 0 && $line[$adjp] < $edgeRpv) {
			$down{$name}{$gene} = 1;
			print DOWN "$tid{$gene}	$gene	$gene{$tid{$gene}}	$entrez{$tid{$gene}}\n";
			if ($entrez{$tid{$gene}} eq '') {
			} else {
				print DOWNE "$entrez{$tid{$gene}}\n";
			}
		}
		if ($logfc{$name}{$gene} > 0 && $line[$adjp] < $edgeRpv) {
			$up{$name}{$gene} = 1;
			print UP "$tid{$gene}	$gene	$gene{$tid{$gene}}	$entrez{$tid{$gene}}\n";
			if ($entrez{$tid{$gene}} eq '') {
			} else {
				print UPE "$entrez{$tid{$gene}}\n";
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
$count = 0;
foreach $i (sort {$a <=> $b} keys %cluster_count) {
	($file) = ($ARGV[1] =~ /(\S+)\.txt/);
	$outfile = $file . ".cluster." . $i . ".txt"; 
	open(OUT, ">$outfile");
	$outfile2 = $file . ".cluster." . $i . ".entrez.txt"; 
	open(OUT2, ">$outfile2");
	foreach $j (@names) {
		$outfile3 = $file . ".cluster." . $i . "." . $j . ".edgeR.down." . $edgeRpv . ".txt";
		$filehandledown{$i}{$j} = "down." . $i . "." . $j;
		open($filehandledown{$i}{$j}, ">$outfile3");
	}
	foreach $j (@names) {
		$outfile4 = $file . ".cluster." . $i . "." . $j . ".edgeR.up." . $edgeRpv . ".txt";
		$filehandleup{$i}{$j} = "up." . $i . "." . $j;
		open($filehandleup{$i}{$j}, ">$outfile4");
	}	
	foreach $j (@names) {
        $outfile5 = $file . ".cluster." . $i . "." . $j . ".edgeR.down.entrez." . $edgeRpv . ".txt";
        $filehandledownentrez{$i}{$j} = "down." . $i . "." . $j . ".entrez";
        open($filehandledownentrez{$i}{$j}, ">$outfile5");
    }       
    foreach $j (@names) { 
        $outfile6 = $file . ".cluster." . $i . "." . $j . ".edgeR.up.entrez." . $edgeRpv . ".txt";
        $filehandleupentrez{$i}{$j} = "up." . $i . "." . $j . ".entrez";
        open($filehandleupentrez{$i}{$j}, ">$outfile6");
    }
	foreach $j (1 .. $cluster_count{$i}) {
		print OUT "$tids[$count]	$gid{$tids[$count]}	$gene{$tids[$count]}	$entrez{$tids[$count]}\n";
		foreach $k (@names) {
			if ($down{$k}{$gid{$tids[$count]}} == 1) {
				print {$filehandledown{$i}{$k}} "$tids[$count]        $gid{$tids[$count]}     $gene{$tids[$count]}    $entrez{$tids[$count]}\n";
				unless ($entrez{$tids[$count]} eq '') {
					print {$filehandledownentrez{$i}{$k}} "$entrez{$tids[$count]}\n";
				}
			}
			if ($up{$k}{$gid{$tids[$count]}} == 1) {
                print {$filehandleup{$i}{$k}} "$tids[$count]        $gid{$tids[$count]}     $gene{$tids[$count]}    $entrez{$tids[$count]}\n";
				unless ($entrez{$tids[$count]} eq '') {
                    print {$filehandleupentrez{$i}{$k}} "$entrez{$tids[$count]}\n";
                }
            }
		}
		unless ($entrez{$tids[$count]} eq '') {
			print OUT2 "$entrez{$tids[$count]}\n";
		}
		$count++;
	}
	close(OUT);
}

		
