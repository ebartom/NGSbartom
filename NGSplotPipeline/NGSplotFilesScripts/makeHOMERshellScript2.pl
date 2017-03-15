$homerShellScript2 = $ARGV[1];
open(OUT, ">$homerShellScript2");
$assembly = $ARGV[2];
opendir(DIR, $ARGV[0]);
@files = grep /entrez/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	if ($filecount == 1) {
		$goDir = $ARGV[0] . "homerGO";
		print OUT "module use /projects/b1025/tools/Modules\n";
		print OUT "module load homer\n";
		print OUT "which findGO.pl\n";
		print OUT "mkdir $goDir\n";
	}	
	($filename) = ($file =~ /(\S+)\.txt/);
	print OUT "# $filename\n";
	$goFileDir = $goDir . "/" . $filename . ".GO";
	print OUT "mkdir $goFileDir\n";
	if ($assembly eq 'hg18' || $assembly eq 'hg19') { 
		$organism = "human";
	} elsif ($assembly eq 'mm9' || $assembly eq 'mm10') {
                $organism = "mouse";
    } elsif ($assembly eq 'dm3' || $assembly eq 'dm6') {
                $organism = "fly";
    } else {
		$organsim = "yeast";
	}
	print OUT "perl /projects/b1025/tools/homer/bin/findGO.pl $ARGV[0]\/$file $organism $goFileDir\n";  
}
