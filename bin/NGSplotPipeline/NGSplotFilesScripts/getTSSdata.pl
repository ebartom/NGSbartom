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
	$chr = $line[1];
	$chr_hash{$chr} = 1;
	$start = $line[2];
	$stop = $line[3];
	$tid = $line[6];
	$gid = $line[4];
	$strand = $line[7];
	$gene = $line[5];
	if ($strand eq '+') {
		$tss = $start;
	} else {
		$tss = $stop;
	}
	push @{$tss{$chr}}, $tss;
	$tid{$chr}{$tss} = $tid; 
	$gid{$tid} = $gid;
	$gene{$tid} = $gene;
	$strand{$tid} = $strand;
}
close(IN);
foreach $i (sort {$a cmp $b} keys %chr_hash) {
	@{$sorted_tss{$i}} = sort {$a <=> $b} @{$tss{$i}}
}
open(IN, $ARGV[1]);
($basename) = ($ARGV[1] =~ /(\S+)\.bed/);
$outfile = $basename . ".nearestGeneList.txt";
open(OUT, ">$outfile");
while($line = <IN>){;
	chomp $line;
	@line = split("\t", $line);
	$chr = $line[0];
	$start = $line[1];
	$stop = $line[2];
	$k27actss{$chr}{$start} = 1;
	$center = int(($stop - $start) / 2) + $start;
   	$max = 300000000;
    	$closest_gene = "";
    	$closest_gene2 = "";
    	foreach $i (@{$sorted_tss{$chr}}) {
    		$distance = abs($center - $i);
    		if ($distance < $max) {
    			$max = $distance;
    			$closest_gene = $gid{$tid{$chr}{$i}};
    			$closest_gene2 = $gene{$tid{$chr}{$i}};
			$closest_gene3 = $tid{$chr}{$i};
    		}
    	}
	$k = $closest_gene;
	$k2 = $closest_gene2;
	$k3 = $closest_gene3;
	$uid = $chr . "_" . $center;
    	print OUT "$uid\t$k\t$k2\t$k3\n";
	$k27acnearestGene{$chr}{$start} = $k;
	$k27acnearestTranscript{$chr}{$start} = $k3;

}
close(IN);
close(OUT);
open(IN, $ARGV[2]);
($basename) = ($ARGV[2] =~ /(\S+)\.bed/);
$outfile = $basename . ".nearestGeneList.txt";
open(OUT, ">$outfile");
while($line = <IN>){;
	chomp $line;
	@line = split("\t", $line);
	$chr = $line[0];
	$start = $line[1];
	$stop = $line[2];
	$k4m3tss{$chr}{$start} = 1;
	$center = int(($stop - $start) / 2) + $start;
   	$max = 300000000;
    	$closest_gene = "";
    	$closest_gene2 = "";
    	foreach $i (@{$sorted_tss{$chr}}) {
    		$distance = abs($center - $i);
    		if ($distance < $max) {
    			$max = $distance;
    			$closest_gene = $gid{$tid{$chr}{$i}};
    			$closest_gene2 = $gene{$tid{$chr}{$i}};
			$closest_gene3 = $tid{$chr}{$i};
    		}
    	}
	$k = $closest_gene;
	$k2 = $closest_gene2;
	$k3 = $closest_gene3;
	$uid = $chr . "_" . $center;
    	print OUT "$uid\t$k\t$k2\t$k3\n";
	$k4m3nearestGene{$chr}{$start} = $k;
	$k4m3nearestTranscript{$chr}{$start} = $k3;

}
close(OUT);
close(IN);
opendir (DIR, $ARGV[3]);
@files = grep /txt/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.edgeR\.txt/);
	$filename_hash{$filename} = 1;
	$infile = $ARGV[3] . $file;
	open(IN, $infile);
	$linecount = 0;
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
		$linecount++;
		$line =~ s/"//g;
		chomp $line;
		@line = split("\t", $line);
		if ($linecount == 1) {
			foreach $i (0 .. $#line) {
				if ($line[$i] =~ /logFC/) {
					$logfc = $i;
				}
			}
			next;
		}
		$gene = $line[0];
		$genelogfc{$filename}{$gene} = $line[$logfc];
	}
	close(IN);
	print "$filename	$linecount\n";
}
close(IN);
opendir (DIR, $ARGV[4]);
@files = grep /edgeR/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.edgeR\.txt/);
	$upbed = $filename . ".k27ac.tss.up.bed";
	$dnbed = $filename . ".k27ac.tss.dn.bed";
	$ncbed = $filename . ".k27ac.tss.nc.bed";
	$upgl = $filename . ".k27ac.tss.up.geneList.txt";
	$dngl = $filename . ".k27ac.tss.dn.geneList.txt";
	$ncgl = $filename . ".k27ac.tss.nc.geneList.txt";
	$upbp = $filename . ".k27ac.tss.up.boxplot.txt";
	$dnbp = $filename . ".k27ac.tss.dn.boxplot.txt";
	$ncbp = $filename . ".k27ac.tss.nc.boxplot.txt";
	open(UPBED, ">$upbed");
	open(DNBED, ">$dnbed");
	open(NCBED, ">$ncbed");
	open(UPGL, ">$upgl");
	open(DNGL, ">$dngl");
	open(NCGL, ">$ncgl");
	open(UPBP, ">$upbp");
	open(DNBP, ">$dnbp");
	open(NCBP, ">$ncbp");
	foreach $i (sort {$a cmp $b} keys %filename_hash) {
		print UPBP "\t$i";
		print DNBP "\t$i";
		print NCBP "\t$i";
	}
	print UPBP "\n";
	print DNBP "\n";
	print NCBP "\n";	
	$infile = $ARGV[4] . $file;
	open(IN, $infile);
	$linecount = 0;
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
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
		$pos = $line[0];
		@pos = split(":", $pos);
		$chr = $pos[0];
		$range = $pos[1];
		@range = split("-", $pos[1]);
		$start = $range[0];
		$end = $range[1];
		unless ($k27actss{$chr}{$start} == 1) {
			next;
		}
		$gene = $k27acnearestGene{$chr}{$start};
		$transcript = $k27acnearestTranscript{$chr}{$start};
		$peaklogfc = $line[$logfc];
		$peakadjp = $line[$adjp];
		if ($peaklogfc > 0 && $peakadjp < 0.05) {
			$k27acuptranscriptcount{$filename}{$transcript}++;
			$k27acupgenecount{$filename}{$gene}++;
			unless ($k27acuptranscriptcount{$filename}{$transcript} > 1) {
				print UPGL "$transcript\n";
			}
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				$genetemp = $gene . "-" . $k27acupgenecount{$filename}{$gene};
				print UPBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print UPBP "\t$genelogfc{$i}{$gene}";
				}
				print UPBP "\n";
			}
			print UPBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		if ($peaklogfc < 0 && $peakadjp < 0.05) {
			$k27acdntranscriptcount{$filename}{$transcript}++;
			$k27acdngenecount{$filename}{$gene}++;
			unless ($k27acdntranscriptcount{$filename}{$transcript} > 1) {
				print DNGL "$transcript\n";
			}
			$genetemp = $gene . "-" . $k27acdngenecount{$filename}{$gene};
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				print DNBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print DNBP "\t$genelogfc{$i}{$gene}";
				}
				print DNBP "\n";
			}
			print DNBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		if ($peakadjp > 0.05) {
			$k27acnctranscriptcount{$filename}{$transcript}++;
			$k27acncgenecount{$filename}{$gene}++;
			unless ($k27acnctranscriptcount{$filename}{$transcript} > 1) {
				print NCGL "$transcript\n";
			}
			$genetemp = $gene . "-" . $k27acncgenecount{$filename}{$gene};
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				print NCBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print NCBP "\t$genelogfc{$i}{$gene}";
				}
				print NCBP "\n";
			}
			print NCBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		
	}
	close(IN);
	close(UPBED);
	close(DNBED);
	close(NCBED);
	close(UPGL);
	close(DNGL);
	close(NCGL);
	close(UPBP);
	close(DNBP);
	close(NCBP);
}
opendir (DIR, $ARGV[5]);
@files = grep /edgeR/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	print "$file\n";
	($filename) = ($file =~ /(\S+)\.edgeR\.txt/);
	$upbed = $filename . ".k4m3.tss.up.bed";
	$dnbed = $filename . ".k4m3.tss.dn.bed";
	$ncbed = $filename . ".k4m3.tss.nc.bed";
	$upgl = $filename . ".k4m3.tss.up.geneList.txt";
	$dngl = $filename . ".k4m3.tss.dn.geneList.txt";
	$ncgl = $filename . ".k4m3.tss.nc.geneList.txt";
	$upbp = $filename . ".k4m3.tss.up.boxplot.txt";
	$dnbp = $filename . ".k4m3.tss.dn.boxplot.txt";
	$ncbp = $filename . ".k4m3.tss.nc.boxplot.txt";
	open(UPBED, ">$upbed");
	open(DNBED, ">$dnbed");
	open(NCBED, ">$ncbed");
	open(UPGL, ">$upgl");
	open(DNGL, ">$dngl");
	open(NCGL, ">$ncgl");
	open(UPBP, ">$upbp");
	open(DNBP, ">$dnbp");
	open(NCBP, ">$ncbp");
	foreach $i (sort {$a cmp $b} keys %filename_hash) {
		print UPBP "\t$i";
		print DNBP "\t$i";
		print NCBP "\t$i";
	}
	print UPBP "\n";
	print DNBP "\n";
	print NCBP "\n";	
	$infile = $ARGV[5] . $file;
	open(IN, $infile);
	$linecount = 0;
	while ($line = <IN>) {
		chomp $line;
		@line = split("\t", $line);
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
		$pos = $line[0];
		@pos = split(":", $pos);
		$chr = $pos[0];
		$range = $pos[1];
		@range = split("-", $pos[1]);
		$start = $range[0];
		$end = $range[1];
		unless ($k4m3tss{$chr}{$start} == 1) {
			next;
		}
		$gene = $k4m3nearestGene{$chr}{$start};
		$transcript = $k4m3nearestTranscript{$chr}{$start};
		$peaklogfc = $line[$logfc];
		$peakadjp = $line[$adjp];
		if ($peaklogfc > 0 && $peakadjp < 0.05) {
			$k4m3uptranscriptcount{$filename}{$transcript}++;
			$k4m3upgenecount{$filename}{$gene}++;
			unless ($k4m3uptranscriptcount{$filename}{$transcript} > 1) {
				print UPGL "$transcript\n";
			}
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				$genetemp = $gene . "-" . $k4m3upgenecount{$filename}{$gene};
				print UPBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print UPBP "\t$genelogfc{$i}{$gene}";
				}
				print UPBP "\n";
			}
			print UPBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		if ($peaklogfc < 0 && $peakadjp < 0.05) {
			$k4m3dntranscriptcount{$filename}{$transcript}++;
			$k4m3dngenecount{$filename}{$gene}++;
			unless ($k4m3dntranscriptcount{$filename}{$transcript} > 1) {
				print DNGL "$transcript\n";
			}
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				$genetemp = $gene . "-" . $k4m3dngenecount{$filename}{$gene};
				print DNBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print DNBP "\t$genelogfc{$i}{$gene}";
				}
				print DNBP "\n";
			}
			print DNBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		if ($peakadjp > 0.05) {
			$k4m3nctranscriptcount{$filename}{$transcript}++;
			$k4m3ncgenecount{$filename}{$gene}++;
			unless ($k4m3nctranscriptcount{$filename}{$transcript} > 1) {
				print NCGL "$transcript\n";
			}
			$switch = 1;
			foreach $i (sort {$a cmp $b} keys %filename_hash) {
				if ($genelogfc{$i}{$gene} eq '') {
					$switch = 0;
				}
			}
			if ($switch == 1) {
				$genetemp = $gene . "-" . $k4m3ncgenecount{$filename}{$gene};
				print NCBP "$genetemp";
				foreach $i (sort {$a cmp $b} keys %filename_hash) {
					print NCBP "\t$genelogfc{$i}{$gene}";
				}
				print NCBP "\n";
			}
			print NCBED "$chr\t$start\t$end\t$pos\t$gene\t+\n";
		}
		
	}
	close(IN);
	close(UPBED);
	close(DNBED);
	close(NCBED);
	close(UPGL);
	close(DNGL);
	close(NCGL);
	close(UPBP);
	close(DNBP);
	close(NCBP);
}