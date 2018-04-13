$homerShellScript = $ARGV[1];
open(OUT, ">$homerShellScript");
$assembly = $ARGV[2];
$promoters = $ARGV[3];
$cgis = $ARGV[4];
opendir(DIR, $ARGV[0]);
@files = grep /bed/, readdir DIR;
close DIR;
$filecount = 0;
foreach $file (@files) {
	$filecount++;
	if ($filecount == 1) {
		$annoDir = $ARGV[0] . "homerAnnotation";
		$motifDir = $ARGV[0] . "homerMotifAnalysis";
		print OUT "module use /projects/b1025/tools/Modules\n";
		print OUT "module load homer\n";
		print OUT "which annotatePeaks.pl\n";
		print OUT "which findMotifsGenome.pl\n";
		print OUT "mkdir $annoDir\n";
		print OUT "mkdir $motifDir\n";
	}	
	($filename) = ($file =~ /(\S+)\.bed/);
	print OUT "# $filename\n";
	$annofile = $annoDir . "/" . $filename . ".anno.txt";
	$annostatsfile = $annoDir . "/" . $filename . ".annoStats.txt";
	unless ($cgis eq '') {
		$cgiDir = $motifDir . "/" . $filename. ".cgi.bg.1000";
		print OUT "mkdir $cgiDir\n";
	}
	$promoterDir = $motifDir . "/" . $filename. ".promoter.bg.1000";
	$genomeDir = $motifDir . "/" . $filename. ".genome.bg.1000";
	$genomeOntologyDir = $annoDir . "/" . $filename . ".genomeOntology";
	$geneOntologyDir = $annoDir . "/" . $filename . ".geneOntology";
	print OUT "mkdir $promoterDir\n";
	print OUT "mkdir $genomeDir\n";
	print OUT "mkdir $genomeOntologyDir\n";
	print OUT "mkdir $geneOntologyDir\n";
	print OUT "perl /projects/b1025/tools/homer/bin/annotatePeaks.pl $ARGV[0]\/$file $assembly -annStats $annostatsfile -go $geneOntologyDir -genomeOntology $genomeOntologyDir > $annofile\n";
	print OUT "perl /projects/b1025/tools/homer/bin/findMotifsGenome.pl $ARGV[0]\/$file $assembly $genomeDir -size 1000\n"; 
	print OUT "perl /projects/b1025/tools/homer/bin/findMotifsGenome.pl $ARGV[0]\/$file $assembly $promoterDir -bg $promoters -size 1000\n";
	unless ($cgis eq '') {
		print OUT "perl /projects/b1025/tools/homer/bin/findMotifsGenome.pl $ARGV[0]\/$file $assembly $cgiDir -bg $cgis -size 1000\n";
	} 
}
