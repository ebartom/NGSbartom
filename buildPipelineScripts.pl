#!/software/activeperl/5.16/bin/perl -w
use Getopt::Long qw(GetOptions);
use List::Util qw(max);
use List::MoreUtils qw(uniq);
use strict;
use utf8;
use warnings;
#use Config::Abstract::Ini;

unless (@ARGV) {
    print "\nUsage: buildPipelineScripts.pl\n\n-config\t\t<configFile>\n-o\t\t<outputDirectory>\n-bs\t\t<baseSpaceDirectory>\n-f\t\t<fastqDirectory>\n";
    print "-bam\t\t<bamDirectory>\n-c\t\tcomparisons.csv file>\n-4C\t\t<4C description file>\n-chip\t\t<ChIP Description file>\n";
    print "-p\t\t<numProcessors>\n-m\t\t<multiMap (1 = allow multi map, 0 = not)>\n";
    print "-a\t\t<aligner>\n-g\t\t<assembly/genome>\n-w\t\t<walltime>\n-t\t\t<RNA|chipseq|4C>\n-account\t<accountName>\n";
    print "-node\t\t<nodeName(eg qnode4144)>\n-scientist\t<initials>\n-s3\t\t<amazon s3 path (if non-standard)>\n-buildBcl2fq\t\t<1|0>\n";
    print "-runBcl2fq\t\t<1|0>\n-runTrim\t\t<1|0>\n-buildAlign\t\t<1|0>\n-runAlign\t\t<1|0>\n-makeTracks\t\t<1|0>\n";
    print "-uploadASHtracks\t<1|0>\n-uploadBAM\t\t<1|0>\n-buildEdgeR\t\t<1|0>\n-runEdgeR\t\t<1|0>\n-buildPeakCaller\t<1|0>\n-runPeakCaller\t\t<1|0>\n";
    print "-buildDiffPeaks\t\t<1|0>\n-runDiffPeaks\t\t<1|0>\n-build4C\t\t<1|0>\n-run4C\t\t\t<1|0>\n";
}

# Set up the environment for analysis.
my $NGSbartom="/projects/p20742/";

my $configFile = "";
my $distToTSS = 2000;
my $upstream = 5000;
my $downstream = 5000;
my $trimString = "TRAILING:30 MINLEN:20";
my $tophatReadMismatch = 2;
my $tophatReadEditDist = 2;
my $tophatMultimap = 5;
my $runPairedEnd = 0;
my $outputDirectory = "";
my $baseSpaceDirectory = "";
my $bamDirectory = "";
my $fastqDirectory = "";
my $sampleSheet = "";
my $comparisons = "";
my $type = "";
my $numProcessors = 4;
my $multiMap = 0;
my $aligner = "";
my $stranded = 1;
my $id = ""; 
my $assembly = "";
my $runBcl2fq = 0;
my $runAlign = 0;
my $runEdgeR = 0;
my $makeTracks = 0;
my $uploadASHtracks = 1;
my $uploadPulmtracks = 0;
my $uploadBAM = 0;
my $buildEdgeR = 0;
my $runTrim = 1;
my $buildBcl2fq = 0;
my $buildAlign = 0;
my $buildGenotyping = 0;
my $runGenotyping = 0;
my $build4C = 0;
my $run4C = 0;
my $buildPeakCaller = 0;
my $runPeakCaller = 0;
my $buildDiffPeaks = 0;
my $runDiffPeaks = 0;
my $walltime = "24:00:00";
my $account = "b1042";
my $queue = "genomics";
my $node = "";
my $scientist = "XXX";
my $s3path = "";
my $fourCdescription = "";
my $chipDescription = "";
my $htseq = 1;
my $bedtools = 0;
my $ngsplot = 0;
my $granges = 0;
my $rsem = 0;
my $genomeBAM = 1;
my $rgString = "";

# S3path is not working properly right now for non-standard s3path values.  This needs fixing.

# Read in the command line arguments.
GetOptions('samplesheet|ss=s' => \$sampleSheet,
	   'config|ini=s' => \$configFile,
	   'outputDirectory|o=s' => \$outputDirectory,
	   'baseSpaceDirectory|bs=s' => \$baseSpaceDirectory,
	   'fastqDirectory|f=s' => \$fastqDirectory,
	   'bamDirectory|bam=s' => \$bamDirectory,
	   'comparisons|c=s' => \$comparisons,
	   '4Cfile|4C=s' => \$fourCdescription,
	   'ChIPfile|chip=s' => \$chipDescription,
	   'type|t=s' => \$type,
	   'aligner|a=s' => \$aligner,
	   'trimString|ts=s' => \$trimString,
	   'assembly|g=s' => \$assembly,
	   'tango|id=s' => \$id,
	   'scientist|s=s' => \$scientist,
	   'stranded|str=i' => \$stranded,
	   'processors|p=i' => \$numProcessors,
	   'walltime|w=s' => \$walltime,
	   'account|acc=s' => \$account,
	   'queue|q=s' => \$queue,
	   'node|n=s' => \$node,
	   's3path|s3=s' => \$s3path,
	   'rgstring|rgs=s' => \$rgString,
	   'tophatMultiMap|tm=i' => \$tophatMultimap,
	   'tophatReadEditDist|ted=i' => \$tophatReadEditDist,
	   'tophatReadMismatch|trm=i' => \$tophatReadMismatch,
	   'multiMap|m=i' => \$multiMap,
	   'buildAlign|ba=i' => \$buildAlign,
	   'buildBcl2fq|bb=i' => \$buildBcl2fq,
	   'runAlign|ra=i' => \$runAlign,
	   'makeTracks|mt=i' => \$makeTracks,
	   'uploadASHtracks|ut=i' => \$uploadASHtracks,
	   'uploadPulmtracks|upt=i' => \$uploadPulmtracks,
	   'uploadBAM|ub=i' => \$uploadBAM,
	   'runTrim|rt=i' => \$runTrim,
	   'runPairedEnd|rpe=i' => \$runPairedEnd,
	   'buildEdgeR|be=i' => \$buildEdgeR,
	   'runEdgeR|re=i' => \$runEdgeR,
	   'buildPeakCaller|bp=i' => \$buildPeakCaller,
	   'runPeakCaller|rp=i' => \$runPeakCaller,
	   'buildDiffPeaks|bdp=i' => \$buildDiffPeaks,
	   'runDiffPeaks|rdp=i' => \$runDiffPeaks,
	   'buildGenotyping|bg=i' => \$buildGenotyping,
	   'runGenotyping|rg=i' => \$runGenotyping,
	   'build4C|b4=i' => \$build4C,
	   'run4C|r4=i' => \$run4C,
	   'runBcl2fq|rb=i' => \$runBcl2fq,
	   'htseq|h=i' => \$htseq,
	   'bedtools|bt=i' => \$bedtools,
	   'ngsplot|np=i' => \$ngsplot,
	   'granges|ash=i' => \$granges,
	   'rsem=i' => \$rsem,
	   'genomeBAM|gb=i' => \$genomeBAM
    ) ;

# if ($configFile){
#     my $abstract = new Config::Abstract::Ini($configFile) or die $!;
#     print STDERR "Contents of Config file:\n$abstract\n";
#     my %PARAMETERS = $abstract->get_entry('PARAMETERS');
#     if (exists($PARAMETERS{'distToTSS'})){ $distToTSS = $PARAMETERS{'distToTSS'};}
#     if (exists($PARAMETERS{'upstream'})){$upstream = $PARAMETERS{'upstream'};}
#     if (exists($PARAMETERS{'downstream'})){$downstream = $PARAMETERS{'downstream'};}
#     if (exists($PARAMETERS{'trimString'})){$trimString = $PARAMETERS{'trimString'};}
#     if (exists($PARAMETERS{'tophatReadMismatch'})){$tophatReadMismatch = $PARAMETERS{'tophatReadMismatch'};}
#     if (exists($PARAMETERS{'tophatReadEditDist'})){$tophatReadEditDist = $PARAMETERS{'tophatReadEditDist'};}
#     if (exists($PARAMETERS{'tophatMultimap'})){$tophatMultimap = $PARAMETERS{'tophatMultimap'};}
#     if (exists($PARAMETERS{'outputDirectory'})){$outputDirectory = $PARAMETERS{'outputDirectory'};}
#     if (exists($PARAMETERS{'baseSpaceDirectory'})){$baseSpaceDirectory = $PARAMETERS{'baseSpaceDirectory'};}
#     if (exists($PARAMETERS{'bamDirectory'})){$bamDirectory = $PARAMETERS{'bamDirectory'};}
#     if (exists($PARAMETERS{'fastqDirectory'})){$fastqDirectory = $PARAMETERS{'fastqDirectory'};}
#     if (exists($PARAMETERS{'sampleSheet'})){$sampleSheet = $PARAMETERS{'sampleSheet'};}
#     if (exists($PARAMETERS{'comparisons'})){$comparisons = $PARAMETERS{'comparisons'};}
#     if (exists($PARAMETERS{'type'})){$type = $PARAMETERS{'type'};}
#     if (exists($PARAMETERS{'numProcessors'})){$numProcessors = $PARAMETERS{'numProcessors'};}
#     if (exists($PARAMETERS{'stranded'})){$stranded = $PARAMETERS{'stranded'};}
    
#     #my $multiMap = 0;
#     #my $aligner = "";
#     #my $id = ""; 
#     #my $assembly = "";
#     #my $runBcl2fq = 0;
#     #my $runAlign = 0;
#     #my $runEdgeR = 0;
#     #my $makeTracks = 0;
#     #my $uploadASHtracks = 1;
#     #my $buildEdgeR = 0;
#     #my $runTrim = 1;
#     #my $buildBcl2fq = 0;
#     #my $buildAlign = 0;
#     #my $build4C = 0;
#     #my $run4C = 0;
#     #my $buildPeakCaller = 0;
#     #my $runPeakCaller = 0;
#     #my $buildDiffPeaks = 0;
#     #my $runDiffPeaks = 0;
#     #my $walltime = "24:00:00";
#     #my $account = "b1025";
#     #my $node = "";
#     #my $scientist = "XXX";
#     #my $s3path = "";
#     #my $fourCdescription = "";
#     #my $chipDescription = "";
# }

if (($uploadPulmtracks == 1) && ($uploadASHtracks == 1)){
    print STDERR "ERROR: Set to upload to both Shilatifard S3 account and Pulmonology S3 account.  Since Shilatifard is default and Pulmonology is manually set, assuming only Pulmonology is desired and resetting flags accordingly.\n";
    $uploadASHtracks = 0;
}
if (($s3path eq "") && ($uploadASHtracks == 1)){ $s3path = "ash-tracks/TANGO/$scientist";}
if (($s3path eq "") && ($uploadPulmtracks == 1)){ $s3path = "m-328-data/$scientist";}
if ($s3path eq ""){ $s3path = "ash-tracks/TANGO/$scientist";}


# Kill the script if output directory, base space directory and type are not specified.
if (($outputDirectory eq "")|| ($type eq "")){
    die "\nERR:  Output Directory and Type must be specified. If Output Directory cannot be found, it will be created.\nPaths should be absolute, NOT relative.\n\n";
}
# If samplesheet is not specified, assume it's inside base space directory 
# with standard name.
if (($sampleSheet eq "") && ($baseSpaceDirectory ne "")){
    $sampleSheet = "$baseSpaceDirectory/SampleSheet.csv";
}
# If aligner is not specified, assume tophat for type RNA and bowtie for type DNA.
if ($aligner eq ""){
    if (($type eq "exome" ) || ($type eq "WGS")){ $aligner = "bwa";}
    if ($type eq "DNA"){ $aligner = "bowtie";}
    if ($type eq "RNA"){ $aligner = "tophat";}
    if ($type eq "chipseq"){ $aligner = "bowtie";}
    if ($type eq "4C"){ $aligner = "bowtie";}
}

# You can build scripts without running them, but not the reverse.
if (($runBcl2fq == 1) && ($buildBcl2fq == 0)){
    $buildBcl2fq = 1;
}
if (($runAlign == 1) && ($buildAlign == 0)){
    $buildAlign = 1;
}
if (($runEdgeR == 1) && ($buildEdgeR == 0)){
    $buildEdgeR = 1;
}
if (($runPeakCaller == 1) && ($buildPeakCaller == 0)){
    $buildPeakCaller = 1;
}
if (($runDiffPeaks == 1) && ($buildDiffPeaks == 0)){
    $buildDiffPeaks = 1;
}
if (($run4C == 1) && ($build4C == 0)){
    $build4C = 1;
}
if (($build4C ==1) && ($fourCdescription eq "")){
    die "ERR:  Cannot de-multiplex 4C samples without file listing primers/viewpoints and indices.\n";
}
if (($buildPeakCaller == 1) && ($chipDescription eq "")){
    die "ERR: Please specify a chipseq description file that pairs ip and input files and specifies narrow or broad peaks.\n";
}
if ($multiMap > 1){
    die "ERR: The MultiMap parameter (-m) is now being used to specify whether the RNAseq tracks should include multi mapped reads.  It should be either 1 or 0.  It used to specify how many alignments tophat should report back for each read.  If this change is causing you problems or consternation, email me at ebartom\@northwestern.edu\n";
}

# Assume that if you want to build Alignment scripts, you want to make tracks too.  This is for backwards compatibility.
if (($buildAlign ==1 ) && ($makeTracks == 0)){
    $makeTracks = 1;
    print STDERR "buildAlign == 1, makeTracks == 0.  Assuming this is an error and switching to makeTracks = 1\n";
}
if ($buildPeakCaller == 1){$makeTracks = 1;}
if (($fastqDirectory ne "") && ($buildBcl2fq == 1)){
    die "ERR:  Cannot de-multiplex from fastq files!  Restart analysis with runBcl2fq and buildBcl2fq = 0\n";
}
my $startFromBAM = 0;
if ($bamDirectory ne ""){ $startFromBAM = 1;}
if ((($bamDirectory ne "") || ($fastqDirectory ne "")) && ($baseSpaceDirectory ne "")) {
    die "ERR:  Please specify either a directory of processed data (fastqs or bams) or a base space directory, not both.\n";
}
if (($bamDirectory ne "") && ($fastqDirectory ne "")) {
    die "ERR:  Please specify either a directory of fastqs or a directory of bams as the starting point of your analysis, not both.\n";
}
if (($bamDirectory ne "") && ($buildAlign == 1)){
    die "ERR:  Cannot align from BAM files!  Restart analysis with runAlign and buildAlign = 0\n";
}
if (($type eq "RNA") && ($granges == 1) && ($bamDirectory eq "")) {
    print STDERR "Cannot run granges without bam file, so running genomic alignment.\n";
    $genomeBAM = 1;
}



# Define a header for the shell scripts.
my $header = "#!/bin/bash\n";
#$header .= "#MSUB -l nodes=1:ppn=$numProcessors\n"; This is now specified within specific scripts, so bcl2fq can have 8 ppn
$header .= "#MSUB -A $account\n";
# Currently I think it only makes sense to specify queue if you are planning on either genomics or genomicsburst, in which case, account should be b1042.
# If you disagree, let me know at ebartom@northwestern.edu 
if ($account eq "b1042"){
    $header .= "#MSUB -q $queue\n";
}
$header .= "#MSUB -l walltime=$walltime\n";
$header .= "#MSUB -m a\n"; # only email user if job aborts
#$header .= "#MSUB -m abe\n"; # email user if job aborts (a), begins (b) or ends (e)
$header .= "#MSUB -j oe\n";
$header .= "#MOAB -W umask=0113\n";
if ($node ne ""){
    $header .= "#MSUB -l hostlist=$node\n";
}
print STDERR "=====================================\n";
print STDERR "Analysis should be initiated either with a base space directory, or a fastq directory, or a bam directory.  Below you can see what was specified in this run.\n";
print STDERR "BaseSpaceDirectory: $baseSpaceDirectory\n";
print STDERR "Fastq Directory: $fastqDirectory\n";
print STDERR "Bam Directory: $bamDirectory\n";
print STDERR "Output Directory: $outputDirectory\n";
print STDERR "Type of Analysis: $type\n";
print STDERR "=====================================\nAnalysis plan:\n";
if ($buildBcl2fq){ print STDERR "Will build scripts for de-multiplexing with Bcl2fq.\n";}
if ($runBcl2fq){ print STDERR "Will run scripts for de-multiplexing with Bcl2fq.\n";}
if ($runAlign && $runTrim) { print STDERR "Will trim fastq files according to trailing quality scores with Trimmomatic ($trimString).\n";}
if ($buildAlign){ print STDERR "Will build scripts for aligning fastq reads with $aligner.\n";}
if ($runAlign){ print STDERR "Will run scripts for aligning fastq reads with $aligner.\n";}
if ($buildPeakCaller){ print STDERR "Will build scripts for calling peaks according to the experimental plan in $chipDescription.\n";}
if ($runPeakCaller){ print STDERR "Will run scripts for calling peaks.\n";}
if ($makeTracks){ print STDERR "Will make tracks appropriate for the UCSC genome browser.\n";
		  if ($uploadASHtracks){ print STDERR "Will upload tracks to Shilatifard account on Amazon S3 (if user has correct credentials).\n";}
		  if ($uploadPulmtracks){ print STDERR "Will upload tracks to Pulmonology account on Amazon S3 (if user has correct credentials).\n";}
		  if ($uploadBAM){ print STDERR "Will upload BAM files to Shilatifard account on Amazon S3 (if user has correct credentials).\n";}
}
if ($buildEdgeR){ print STDERR "Will build scripts for finding differentially expressed genes from RNAseq data.\n";}
if ($runEdgeR){ print STDERR "Will run scripts for finding differentially expressed genes from RNAseq data.\n";}		  
if ($buildDiffPeaks) { print STDERR "Will build scripts to find differentially expressed peaks based on $chipDescription and $comparisons.\n";}
if ($runDiffPeaks) { print STDERR "Will run scripts to find differentially expressed peaks.\n";}
if ($build4C){ print STDERR "Will build scripts for analyzing 4C data, starting with mock-demultiplexed fastq files.\n";}
if ($run4C){ print STDERR "Will run scripts for analyzing 4C data.\n";}
if ($buildGenotyping){ print STDERR "Will build scripts for genotyping samples; currently only implemented for RNA.\n";}
if ($runGenotyping){ print STDERR "Will run scripts for genotyping samples.\n";}
print STDERR "=====================================\n";

# Define references.
my (%bowtieIndex,%txIndex,%txdbfile,%bwaIndex,%gff,%exonbed,%rsemTx);
my (%gatkRef,%knownSNPsites,%knownIndelsites);

$bowtieIndex{"hg38"} = "$NGSbartom/anno/bowtie_indexes/hg38";
$bwaIndex{"hg38"} = "$NGSbartom/anno/bwa_indexes/hg38.fa";
$txIndex{"hg38"} ="$NGSbartom/anno/tophat_tx/hg38.Ens_78.remap";
$txdbfile{"hg38"} = "$NGSbartom/anno/Txdb/hsapiens_gene_ensembl_Ens78.txdb";
$exonbed{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.exons.bed";
$gff{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.cuff.gtf";
$rsemTx{"hg38"} = "$NGSbartom/anno/rsemTx/hg38.Ens_78";
$gatkRef{"hg38"} = "$NGSbartom/anno/picardDict/hg38.fa";
$knownSNPsites{"hg38"} = "$NGSbartom/anno/picardDict/1000G_phase1.snps.high_confidence.hg38.vcf";
$knownIndelsites{"hg38"} = "$NGSbartom/anno/picardDict/Mills_and_1000G_gold_standard.indels.hg38.vcf";

$bowtieIndex{"hg19"} = "$NGSbartom/anno/bowtie_indexes/hg19";
$bwaIndex{"hg19"} = "$NGSbartom/anno/bwa_indexes/hg19.fa";
$txIndex{"hg19"} ="$NGSbartom/anno/tophat_tx/hg19.Ens_72.remap";
$txdbfile{"hg19"} = "$NGSbartom/anno/Txdb/hsapiens_gene_ensembl_Ens72.txdb";
$exonbed{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_72/hg19.Ens_72.exons.bed";
$gff{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_72/hg19.Ens_72.cuff.gtf";
$rsemTx{"hg19"} = "$NGSbartom/anno/rsemTx/hg19.Ens_72";
$gatkRef{"hg19"} = "$NGSbartom/anno/picardDict/hg19.fa";
$knownSNPsites{"hg19"} = "$NGSbartom/anno/picardDict/1000G_phase1.snps.high_confidence.hg19.sites.noContigs.vcf";
$knownIndelsites{"hg19"} = "$NGSbartom/anno/picardDict/1000G_phase1.indels.hg19.sites.noContigs.vcf";

$bowtieIndex{"dm3"} = "$NGSbartom/anno/bowtie_indexes/dm3";
$bwaIndex{"dm3"} = "$NGSbartom/anno/bwa_indexes/dm3.fa";
$txIndex{"dm3"} = "$NGSbartom/anno/tophat_tx/dm3.Ens_74.cuff";
$txdbfile{"dm3"} = "$NGSbartom/anno/Txdb/dmelanogaster_gene_ensembl_Ens74.txdb";
$exonbed{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.exons.bed";
$gff{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.cuff.gtf";
$rsemTx{"dm3"} = "$NGSbartom/anno/rsemTx/dm3.Ens_74";

$bowtieIndex{"mm10"} = "$NGSbartom/anno/bowtie_indexes/mm10";
$bwaIndex{"mm10"} = "$NGSbartom/anno/bwa_indexes/mm10.fa";
$txIndex{"mm10"} = "$NGSbartom/anno/tophat_tx/mm10.Ens_78.cuff";
$txdbfile{"mm10"} = "$NGSbartom/anno/Txdb/mmusculus_gene_ensembl_Ens78.txdb";
$exonbed{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.exons.bed";
$gff{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf";
$rsemTx{"mm10"} = "$NGSbartom/anno/rsemTx/mm10.Ens_78";

$bowtieIndex{"mm9"} = "$NGSbartom/anno/bowtie_indexes/mm9";
$bwaIndex{"mm9"} = "$NGSbartom/anno/bwa_indexes/mm9.fa";
$txIndex{"mm9"} = "$NGSbartom/anno/tophat_tx/mm9.Ens_67.remap";
$txdbfile{"mm9"} = "$NGSbartom/anno/Txdb/mmusculus_gene_ensembl_Ens67.txdb";
$exonbed{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.exons.bed";
$gff{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.cuff.gtf";
$rsemTx{"mm9"} = "$NGSbartom/anno/rsemTx/mm9.Ens_67";

$bowtieIndex{"sacCer3"} = "$NGSbartom/anno/bowtie_indexes/sacCer3";
$bwaIndex{"sacCer3"} = "$NGSbartom/anno/bwa_indexes/sacCer3.fa";
$txIndex{"sacCer3"} = "$NGSbartom/anno/tophat_tx/sacCer3.Ens_72.remap";
#$txdbfile{"sacCer3"} = "$NGSbartom/anno/Txdb/scerevisiae_gene_ensembl_Ens72.txdb";
#$txIndex{"sacCer3"} = "$NGSbartom/anno/tophat_tx/sacCer3.Ens_78.remap";
$txdbfile{"sacCer3"} = "$NGSbartom/anno/Txdb/scerevisiae_gene_ensembl_Ens78.txdb";
$exonbed{"sacCer3"} = "$NGSbartom/anno/Ens/sacCer3.Ens_78/sacCer3.Ens_78.exons.bed";
$gff{"sacCer3"} = "$NGSbartom/anno/Ens/sacCer3.Ens_78/sacCer3.Ens_78.cuff.gtf";
$rsemTx{"sacCer3"} = "$NGSbartom/anno/rsemTx/sacCer3.Ens_72";

if (-d "$outputDirectory") {
    print STDERR "Output Directory found.\n";
}else {system("mkdir $outputDirectory");print STDERR "Output directory created.\n";}

my $cmd;
if ($sampleSheet ne ""){
# check if SampleSheet contains windows or mac carriage returns, and remove them, if found.
    my $wrongCarriageReturn = `grep \"\\r\" $sampleSheet`;
#print STDERR "Carriage return test = $wrongCarriageReturn\n";
    if ($wrongCarriageReturn ne ""){
	&datePrint("Sample sheet contained windows or mac carriage returns.  Cleaning it.");
	$cmd = "mv $sampleSheet $sampleSheet.old\nperl -pe \"s\/\\r\\n\/\\n\/g\" $sampleSheet.old | perl -pe  \"s\/\\r\/\\n\/g\" > $sampleSheet.fixed.txt\nmv $sampleSheet.fixed.txt $sampleSheet\n";
	system($cmd); 
    }
}

if ( $fourCdescription ne ""){
# check if SampleSheet contains windows or mac carriage returns, and remove them, if found.
    my $wrongCarriageReturn = `grep \"\\r\" $fourCdescription`;
#print STDERR "Carriage return test = $wrongCarriageReturn\n";
    if ($wrongCarriageReturn ne ""){
	&datePrint("Sample sheet contained windows or mac carriage returns.  Cleaning it.");
	$cmd = "mv $fourCdescription $fourCdescription.old\nperl -pe \"s\/\\r\\n\/\\n\/g\" $fourCdescription.old | perl -pe  \"s\/\\r\/\\n\/g\" > $fourCdescription.fixed.txt\nmv $fourCdescription.fixed.txt $fourCdescription\n";
	system($cmd); 
    }
}

if ($comparisons ne ""){
# check if Comparison File contains windows or mac carriage returns, and remove them, if found.
    my $wrongCarriageReturn = `grep \"\\r\" $comparisons`;
#print STDERR "Carriage return test = $wrongCarriageReturn\n";
    if ($wrongCarriageReturn ne ""){
	&datePrint("Comparison file contained windows or mac carriage returns.  Cleaning it.");
	$cmd = "mv $comparisons $comparisons.old\nperl -pe \"s\/\\r\\n\/\\n\/g\" $comparisons.old | perl -pe  \"s\/\\r\/\\n\/g\" | perl -pe \"s/ //g\" > $comparisons.fixed.txt\nmv $comparisons.fixed.txt $comparisons\n";
	system($cmd); 
    }
}

if ($chipDescription ne ""){
# check if ChIP description File contains windows or mac carriage returns, and remove them, if found.
    my $wrongCarriageReturn = `grep \"\\r\" $chipDescription`;
#print STDERR "Carriage return test = $wrongCarriageReturn\n";
    if ($wrongCarriageReturn ne ""){
	&datePrint("Comparison file contained windows or mac carriage returns.  Cleaning it.");
	$cmd = "mv $chipDescription $chipDescription.old\nperl -pe \"s\/\\r\\n\/\\n\/g\" $chipDescription.old | perl -pe  \"s\/\\r\/\\n\/g\" | perl -pe \"s\/ //g\" > $chipDescription.fixed.txt\nmv $chipDescription.fixed.txt $chipDescription\n";
	system($cmd); 
    }
}

# If buildBcl2fq == 1, then create a shell script for Bcl2fq (if runBcl2fq == 1, then submit the job and wait for it to finish).
if ($buildBcl2fq == 1){
    &datePrint("Creating shell script for Bcl2fq");
    my $shScript = "$baseSpaceDirectory\/runBcl2fq.sh";
    my $bclFqProcessors = max($numProcessors,8);
    open(SH,">$shScript");
    print SH $header;
    print SH "#MSUB -l nodes=1:ppn=$bclFqProcessors\n";
    print SH "#MSUB -N bcl2fastq\n";
    print SH "module load bcl2fastq/2.17.1.14\n";
    print SH "bcl2fastq -R $baseSpaceDirectory -r $numProcessors -d $numProcessors -p $numProcessors -w $numProcessors\n";
    # If runBcl2fq == 1, then run the shell script (only works if buildBcl2fq == 1)
    if ($runBcl2fq == 1){
	&datePrint("Running Bcl2fq job and waiting for it to finish.");
	# Submit the job, saving the job id.
	my $result = `msub $shScript`;
	$result =~ s/\s//g;
	my $jobfinished = "no";
	# Wait until the job is no longer in qstat.
	&datePrint("Waiting for job $result to finish.");
	# Check qstat every 30 seconds, adding a "." every time you check.
	until ($jobfinished eq "Completed"){
	    $jobfinished = `checkjob $result | grep ^State:`;
	    sleep(30);
	    print STDERR ".";
	    if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
	    print STDERR "$jobfinished";
	}
	print STDERR "\n";
	&datePrint("Job $result done.  Continuing.");
    }
    $cmd = "cp $baseSpaceDirectory\/runBcl2fq.sh $outputDirectory\/";
    system($cmd);
}
my %scientists;

my(%fastqs,%reference,%samples,%sampleIDs,$sample_project,$sample_plate);
$sample_plate = $scientist;
if (($sampleSheet ne "")){
# Read in the sample sheet
    &datePrint("Reading in $sampleSheet and printing Sample_Report");
    open (IN,$sampleSheet);
    my $flag="header";
    my($sample_ID,$sample_name,$assembly,$I7_Index_ID,$index,$description,$index2,$I5_Index_ID,$stuff);
    my $sampleNum = 0;
# Create output file for Sample_Report.
    open(OUT,">$outputDirectory/Sample_Report.csv");
# Print labels to output file.
    print OUT "Fastq,Sample_Name,Assembly\n";
    while(<IN>){
	chomp $_;
	if (($flag eq "data") && ($_ !~ /^Sample_ID/)){
	    # Read in the metadata from the sample sheet.
	    ($sample_ID,$sample_name,$sample_plate,$assembly,$I7_Index_ID,$index,$sample_project,$description,$stuff) = split(/\,/,$_);
	    $sampleNum++;
	    if ($type ne "4C"){
		$sample_name =~ s/\_/\-/g;
		$sample_name =~ s/\./\-/g;
	    }
	    $sampleIDs{$sample_name} = $sample_ID;
	    if ($description =~ /^[ACTG]+$/) { # If description is a nucleotide sequence, this sample has two indices.  The project name will be in the stuff variable.
		$sample_project = $stuff;
		&datePrint("This sample has two indices.  Using $stuff as project name.");
	    }
	    # Store the assembly for the sample and project.
	    $reference{$sample_name}=$assembly;
	    $reference{$sample_project} = $assembly;
	    $scientists{$sample_project} = $sample_plate;
	    print "$sample_name\tREF1:$reference{$sample_name}\t$bowtieIndex{$assembly}\n";
	    if ($type ne "4C"){
		my $fastq = "";
		foreach my $lane ("L001","L002","L003","L004"){
		    # Build fastq file names.
		    # NB:  Right now it assumes single end.
		    $fastq = "$sample_name\_S$sampleNum\_$lane\_R1\_001.fastq.gz";
		    print STDERR "Looking for $fastq\n";
		    print STDERR "Looking for $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq and $outputDirectory\/$sample_project\/fastq\/$fastq\n";
		    # Check if fastq exists in new subdirectory or old.
		    if ((-e "$baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq") || (-e "$outputDirectory\/$sample_project\/fastq\/$fastq")){
			#		print STDERR "$fastq exists in $baseSpaceDirectory\n";
		    } else { 
			# If you still can't find it, try all underscores (like in SampleSheet), again in old subdirectory or new
			$fastq = "$sample_name\_S$sampleNum\_$lane\_R1\_001.fastq.gz";
			$fastq =~ s/-/_/g;
			print STDERR "Looking for $fastq\n";
			if ((-e "$baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq") || (-e "$outputDirectory\/$sample_project\/fastq\/$fastq")) {
			    print STDERR "$fastq exists\n";
			} else { 
			    # If you can't find Fastq file anywhere, then die.
			    # This would indicate an error in the path or filename, or that the bcl2fq job failed.  It could also indicate that the index was wrong in the sample sheet, and that no reads were assigned to a specific sample.
			    die "ERR:  Cannot find fastq file!\n";
			}
		    }
		    # Write the sample report file.
		    print OUT "$fastq,$sample_name,$assembly\n";
		    # Create a hash of all fastq files for a given sample_name.
		    if (!exists($fastqs{$sample_name})){
			$fastqs{$sample_name} = "$outputDirectory\/$sample_project\/fastq\/$fastq";
		    }else {$fastqs{$sample_name} .= ",$outputDirectory\/$sample_project\/fastq\/$fastq";}
		    print STDERR "$sample_name\t$fastqs{$sample_name}\n";
		}
	    }
	    # Create a hash of all samples in a give sample project (TANGO/MOLNG)
	    if (!exists($samples{$sample_project})){
		$samples{$sample_project} = $sample_name;
	    }else {$samples{$sample_project} .= ",$sample_name";}
	}
	# Check for [Data] in the SampleSheet to see when the data actually starts.
	if ($_ =~ /^\[Data\]/){
	    $flag = "data";
	}
    }
    close OUT;
} elsif ($sampleSheet eq ""){
    my ($sample_name);
    if (($assembly eq "") || (($fastqDirectory eq "")&& ($startFromBAM == 0))){
	die "ERR:  If SampleSheet is not specified, assembly and fastqDirectory or bamDirectory must be\n";
    } elsif ($fastqDirectory ne ""){
	&datePrint("Looking for Fastq files in $fastqDirectory.");
	my $fastqlist = "";
	#	$fastqlist = system("ls $fastqDirectory\/*.fastq.*gz");
	$fastqlist = `ls $fastqDirectory\/*.fastq.gz`;
#	print STDERR $fastqlist;
	if ($fastqlist eq ""){
	    $fastqlist = `ls $fastqDirectory\/*.fq.gz`;
	}
	if ($fastqlist eq ""){
	    $fastqlist = `ls $fastqDirectory\/*.fastq.tgz`;
	}
	if ($fastqlist eq ""){
	    $fastqlist = `ls $fastqDirectory\/*.fastq`;
	}
	if ($fastqlist eq ""){
	    $fastqlist = `ls $fastqDirectory\/*.fq`;
	}
	my @fastqlist = split(/\s+/,$fastqlist);
	&datePrint("Found @fastqlist");
	my $project_name = "thisProject";
	if ($fastqDirectory =~ /([\w\-\_\.]+)\/fastq\/?$/){
	    $project_name = $1;
	} elsif ($fastqDirectory =~ /([\w\-\_\.]+)\/?$/){
	    $project_name = $1;
	    $project_name =~ s/.fastqs//g;
	    $project_name =~ s/.fastq//g;
	    $project_name =~ s/.seqfiles//g;
	}
	$reference{$project_name}=$assembly;
	&datePrint("Project name is $project_name");
#	print STDERR "$project_name @fastqlist\n";
	foreach my $fastq (@fastqlist){
#	    print STDERR "Fastq: \"$fastq\"\n";
	    if (($fastq =~ /\/?([\w\-\d\_\.]+)\_S\d/) || 
		#		($fastq =~ /\/?([\w\-\d\_\.]+)FastqRd/) ||
		($fastq =~ /\/?([\w\-\d\_\.]+)\_R\d/) ||
		($fastq =~ /\/?([\w\-\d\_\.]+).fastq.t?gz/) 
		){
		$sample_name = $1;
		&datePrint("Sample name is $sample_name");
		if (!exists($fastqs{$sample_name})){
		    $fastqs{$sample_name} = "$fastq";
#		    &datePrint($fastqs{$sample_name});
		    $reference{$sample_name}=$assembly;
		}else {$fastqs{$sample_name} .= ",$fastq";}
		# Create a hash of all samples in a give sample project (TANGO/MOLNG)
		if (!exists($samples{$project_name})){
		    $samples{$project_name} = $sample_name;
		}else {$samples{$project_name} .= ",$sample_name";}
	    } 
	}
    } elsif ($startFromBAM == 1){
	&datePrint("Looking for bam files in $bamDirectory.");
	my $bamlist = "";
#	$bamlist = system("ls $bamDirectory\/*.bam");
	$bamlist = `ls $bamDirectory\/*.bam`;
	print STDERR $bamlist;
	my @bamlist = split(/\n/,$bamlist);
	&datePrint("Found @bamlist");
	my $project_name = "thisProject";
	if ($bamDirectory =~ /\/?([\w\-\.\_]+)\/bam\/?$/){
	    $project_name = $1;
	}elsif ($bamDirectory =~ /\/?([\w\-\.\_]+)\/?$/){
	    $project_name = $1;
	    $project_name =~ s/.bams//g;
	    $project_name =~ s/.bam//g;
	    $project_name =~ s/.bamfiles//g;
	}
	$reference{$project_name}=$assembly;
	&datePrint("Project name is $project_name");
#	print STDERR "$project_name @bamlist\n";
	foreach my $bam (@bamlist){
#	    print STDERR "Bam: \"$bam\"\n";
	    if ($bam =~ /\/([\w\-\d\.]+).bam$/){
		$sample_name = $1;
		&datePrint("Sample name is $sample_name");
		if (!exists($reference{$sample_name})){
		    $reference{$sample_name}=$assembly;
		}
		# Create a hash of all samples in a given sample project (TANGO/MOLNG)
		if (!exists($samples{$project_name})){
		    $samples{$project_name} = $sample_name;
		}else {$samples{$project_name} .= ",$sample_name";}
	    } 
	}
    }
}

 
if (($type eq "4C") && ($build4C == 1)){
    my $maxPrimerMismatch = 1;	
    if ($fourCdescription ne ""){
	&datePrint("4C experiment described in $fourCdescription");
	$cmd = "mkdir $outputDirectory\/$sample_project\n";
	$cmd .= "mkdir $outputDirectory\/$sample_project\/scripts\n";
	$cmd .= "mkdir $outputDirectory\/$sample_project\/fastq\n";
	system($cmd);
	open(FCSH,">$outputDirectory\/$sample_project\/scripts\/run_4C_demultiplex.sh");
	print FCSH "$header";
	print FCSH "#MSUB -N 4Cdemultiplex\n";
	print FCSH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	print FCSH "export PATH=\$PATH:$NGSbartom/tools/\n";
	print FCSH "date\n";
	print FCSH "\n# Copy raw reads into fastq directory and de-compress them.\n";
	print FCSH "cp $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/Undetermined*.fastq.gz $outputDirectory\/$sample_project\/fastq\/\n";
	print FCSH "gunzip $outputDirectory\/$sample_project\/fastq\/Undetermined*fastq.gz\n";
	print FCSH "date\n";
	print FCSH "\n# De-multiplex reads using $fourCdescription table\n";
	print FCSH "perl $NGSbartom/tools/split4CwithTable.pl $outputDirectory\/$sample_project\/fastq\/ $fourCdescription $maxPrimerMismatch\n";
	print FCSH "date\n";
	print FCSH "\n# Clean up extra files.\n";
	print FCSH "ls $outputDirectory\/$sample_project\/fastq\/*\_S*\_L00*\_R*fastq\n";
	print FCSH "rm $outputDirectory\/$sample_project\/fastq\/*\_S*\_L00*\_R*fastq\n";
	print FCSH "date\n";
	close(FCSH);
	if (($run4C ==1)){
	    &datePrint("Submitting job to de-multiplex 4C samples");
	    my $result = `msub $outputDirectory/$sample_project/scripts/run_4C_demultiplex.sh`;
	    $result =~ s/\s+//g;
	    my $jobfinished = "no";
	    # Wait until the job is Complete.
	    &datePrint("Waiting for job $result to finish. (each . = 30 seconds)");
	    # Check qstat every 30 seconds, adding a "." every time you check.
	    until ($jobfinished eq "Completed"){
		$jobfinished = `checkjob $result | grep ^State:`;
		sleep(30);
		print STDERR ".";
		if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
		print STDERR "$jobfinished";
	    }
	    print STDERR "\n";
	    &datePrint("Job $result done.  Continuing.");
	}
	my $sampleList = $samples{$sample_project};
	my @samples = uniq(split(/\,/,$sampleList));
	open(OUT,">$outputDirectory/Sample_Report.csv");
	# Print labels to output file.
	print OUT "Fastq,Sample_Name,Assembly\n";
	foreach my $sample (@samples){
#	    print STDERR "$sample\n";
	    $sample =~ s/\_S\d+\_L\d+\_R\d+\_\d+//g;
	    my $fastq = "$sample\.m$maxPrimerMismatch\.fastq.gz";
#	    print STDERR "$sample\t$fastq\n";
	    if (exists($fastqs{$sample})){
		#skip;
	    } else {
		$fastqs{$sample} = "$outputDirectory\/$sample_project\/fastq\/$fastq";
#		print "$fastqs{$sample}\n";
		if (-e $fastqs{$sample}){
		    print STDERR "$fastq exists\n";
		}
		print OUT "$fastq,$sample,$assembly\n";
	    }
	}
    }
}



# Create the directory structure, and move fastq files to a project specific directory (by TANGO/MOLNG) in the top level of the base space directory.
&datePrint("Setting up directory structure and maybe moving fastq files to project sub-directory within $baseSpaceDirectory");
foreach my $project (keys(%samples)){
    if ($scientists{$project} ne ""){ 
	$scientist = $scientists{$project};
#	if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	$s3path = "ash-tracks/TANGO/$scientist";
    }
    my $cmd = "mkdir $outputDirectory\/$project\n";
    $cmd .= "mkdir $outputDirectory\/$project/scripts\n";    
    $cmd .= "mkdir $outputDirectory\/$project\/bam\n";
    print STDERR "$cmd\n";
    system($cmd);
    my @sampleSet = uniq(split(/\,/,$samples{$project}));
    print STDERR "SampleSet: @sampleSet\n";
    my @fastqSet = split(/\,/,$fastqs{$sampleSet[0]});
    print STDERR "FastqSet for $sampleSet[0]:  @fastqSet\n";
    if (($baseSpaceDirectory ne "") && ($type ne "4C")){
	$cmd = "mkdir $outputDirectory\/$project\n";
	$cmd .= "mkdir $outputDirectory\/$project\/fastq\n";
	system($cmd);
	&datePrint ("Checking for $fastqSet[0] to decide whether $project fastqs should be moved.");
	# If the fastq file does not exist in the place the aligner will expect it, move it there.
	if (!(-e "$fastqSet[0]")){
	    &datePrint("Moving fastqs to $outputDirectory\/$project\/fastq\/");
	    $cmd = "mv $baseSpaceDirectory/Data/Intensities/BaseCalls/$project/*/*.fastq.gz $outputDirectory\/$project\/fastq\/";
	    print STDERR "$cmd\n";
	    system($cmd);
	} else { &datePrint("Not over-writing fastqs in $outputDirectory\/$project\/fastq");}
    }
}

# If the aligner is tophat, create the shell scripts to run tophat on all fastqs, one sample at a time.
if (($buildAlign == 1) && ($aligner eq "tophat")){
    &datePrint("Creating Tophat Alignment shell scripts.");
    my @samples;
    # Foreach project (TANGO/MOLNG):
    foreach my $project (keys(%samples)){
	if ($scientists{$project} ne ""){ 
	    $scientist = $scientists{$project};
	    #if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	    $s3path = "ash-tracks/TANGO/$scientist";
	}
	# Make a directory for the output.
	$cmd = "mkdir $outputDirectory\/$project\/Tophat_aln";
	system($cmd);
	@samples = uniq(split(/\,/,$samples{$project}));
	
	# Foreach sample within the project:
	foreach my $sample (@samples){
	    # Create a shell script to run tophat on all fastqs for the sample at the same time.
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_align.sh";
	    &datePrint("Printing to $shScript");
	    open (SH,">$shScript");
	    print SH "$header";
	    print SH "#MSUB -N $sample\_tophat\n";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    print SH "export PATH=\$PATH:$NGSbartom/tools/\n";
	    print SH "module load bowtie2/2.2.6\n";
	    print SH "module load tophat/2.1.0\n";
	    print SH "module load samtools/1.2\n";
	    print SH "module load boost/1.56.0\n";
	    print SH "module load gcc/4.8.3\n";
#	    print SH "module load boost/1.57.0\n\n";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    if ($runTrim == 1){
		print SH "module load java/jdk1.8.0_25\n";
		print SH "\n# Make Directory for FastQC reports\n";
		print SH "mkdir $outputDirectory\/$project\/fastqc\n";
		print SH "mkdir $outputDirectory\/$project\/fastq\n";
		my @newfastqs = ();
		foreach my $fastq (@fastqs){
		    my $fastqname = "";
		    my $newfastq = "";
		    if ($fastq =~ /\/?([\w\d\-\_\.]+\.fastq\.t?gz$)/){
			$fastqname = $1;
		    } elsif ($fastq =~ /\/?([\w\d\-\_\.]+\.fastq$)/){
			$fastqname = $1;
		    }
		    $newfastq = "$outputDirectory\/$project\/fastq\/$fastqname";
#		    print STDERR "Fastq: $fastq\nNewFastq: $newfastq\nFastqname = $fastqname\n";
		    print SH "\n# Trim poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $outputDirectory\/$project\/fastq\/$fastqname.trimmed $trimString\n";
		    print SH "gzip $outputDirectory\/$project\/fastq\/$fastqname.trimmed\n";
		    if ($newfastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.t?gz$/){
			if (-f "$1\_fastqc.html"){
			    print SH "\# FastQC file already exists\n";
			} else {
			    print SH "# Running FastQC to assess read quality.\n";
			    print SH "date\n$NGSbartom/tools/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $fastq $fastq.trimmed.gz\n";
			}
		    }
		    print SH "mv $newfastq.trimmed.gz $newfastq\n";
		    print SH "date\n\n";
		    push(@newfastqs,$newfastq);

		}
		$fastqs{$sample} = "@newfastqs";
		$fastqs{$sample} =~ s/\s/\,/g;
#		&datePrint("New fastqs for sample $sample are @newfastqs, and $fastqs{$sample}");
	    }
	    if ($runPairedEnd == 0){			   
		if ($genomeBAM == 1){
		   
		    print SH "\n# Run Tophat to align data for single end data.\n";
		    if ($buildGenotyping == 1){
			# If genotyping, add read groups.  These could be made more accurate.
			if ($rgString eq ""){
			    $rgString = "--rg-sample $sample --rg-id $sample --rg-library $sample --rg-description $sample --rg-platform-unit nextseq --rg-center ASH --rg-platform nextseq";
			}
		    } 
		    print SH "\ntophat --no-novel-juncs --read-mismatches $tophatReadMismatch --read-edit-dist $tophatReadEditDist --num-threads $numProcessors --max-multihits $tophatMultimap $rgString --transcriptome-index $txIndex{$reference{$sample}} -o \$TMPDIR\/$sample $bowtieIndex{$reference{$sample}} $fastqs{$sample} >& $outputDirectory\/$project\/bam\/$sample.tophat.log\n";
		}
	    } elsif ($runPairedEnd == 1){
		my @read1fastqs;
		my @read2fastqs;
		my @read3fastqs;
		my @fastqs = split(/\,/,$fastqs{$sample});
		foreach my $fastq (@fastqs){
#		    print STDERR "Looking for Read number in $fastq\n";
		    if (($fastq =~ /\_R1\_?\.?/)  ){
#			|| ($fastq =~ /Rd1/)){
			push (@read1fastqs,$fastq);
		    } elsif (($fastq =~ /\_R2\_?\.?/)){
#			|| ($fastq =~ /Rd2/)){
			push (@read2fastqs,$fastq);
		    } elsif ($fastq =~ /\_R3\_?\.?/){
			push (@read3fastqs,$fastq);
		    } else {
			print STDERR "ERR: Could not find Read Number!\n";
		    }
		}
		@read1fastqs = sort(@read1fastqs);
		@read2fastqs = sort(@read2fastqs);
		@read3fastqs = sort(@read3fastqs);
		my $read1fastqs = "@read1fastqs";
		my $read2fastqs = "@read2fastqs";
		my $read3fastqs = "@read3fastqs";
		if (length($read3fastqs)>length($read2fastqs)){
		    $read2fastqs = $read3fastqs;
		}
		$read1fastqs =~ s/\ /,/g;
		$read2fastqs =~ s/\ /,/g;
		if ($genomeBAM == 1){
		    print SH "\n# Run Tophat to align data for paired end data.\n";
		    print SH "\ntophat --no-novel-juncs --read-mismatches $tophatReadMismatch --read-edit-dist $tophatReadEditDist --num-threads $numProcessors --max-multihits $tophatMultimap --transcriptome-index $txIndex{$reference{$sample}} -o \$TMPDIR\/$sample $bowtieIndex{$reference{$sample}} $read1fastqs $read2fastqs >& $outputDirectory\/$project\/bam\/$sample.tophat.log\n";
		} 
		if ($rsem == 1){ # AND runpaired = 1
		    print SH "\n# Align fastqs to transcriptome for RSEM\n";
		    print SH "module load bowtie/1.1.2 \n";
		    print SH "export PATH=\$PATH:$NGSbartom/tools/RSEM-1.2.30/\n";
		    print SH "date\n";
		    print SH "\n# First prepare fastqs\n";
		    #		print SH "mkdir $outputDirectory/$project/fastq/\n";
		    $read1fastqs =~ s/\,/\ /g;
		    $read2fastqs =~ s/\,/\ /g;
		    my $strandstring = "";
		    #if ($stranded == 1) { $strandstring = "--strand-specific";}
		    if ($stranded == 1) { $strandstring = "--forward-prob 0";}
		    print SH "gunzip $read1fastqs\n";
		    print SH "gunzip $read2fastqs\n";
		    $read1fastqs =~ s/.gz//g;
		    $read2fastqs =~ s/.gz//g;
#		    print SH "echo $read1fastqs\n";
#		    print SH "echo $read2fastqs\n";
		    print SH "cat $read1fastqs > $outputDirectory/$project/fastq/$sample.read1.fastq\n";
		    print SH "cat $read2fastqs > $outputDirectory/$project/fastq/$sample.read2.fastq\n";
		    print SH "gzip $read1fastqs &\n";		    
		    print SH "gzip $read2fastqs &\n";
		    print SH "\n# Then calculate expression for $sample.\n";
		    
		    print SH "rsem-calculate-expression --paired-end $outputDirectory/$project/fastq/$sample.read1.fastq $outputDirectory/$project/fastq/$sample.read2.fastq $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample --no-bam-output -p $numProcessors $strandstring --estimate-rspd >& $outputDirectory/$project/bam/$sample.rsem.log\n";
		    print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n";
		    print SH "ls $outputDirectory/$project/fastq/$sample.read*.fastq\n";
		    print SH "rm $outputDirectory/$project/fastq/$sample.read1.fastq\n";
		    print SH "rm $outputDirectory/$project/fastq/$sample.read2.fastq\n";
		}
	    }
	    print SH "date\n";
	    print SH "\nrsync -av \$TMPDIR\/$sample/* $outputDirectory\/$project\/Tophat_aln\/$sample/\n";
#	    print SH "\nmv \$TMPDIR\/$project\_$sample\_tophatOut $outputDirectory\/$project\/Tophat_aln\/$sample\n";
	    print SH "date\n";
	    print SH "ln -s $outputDirectory\/$project\/Tophat_aln\/$sample\/accepted_hits.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n";
	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    if (($rsem == 1)&& ($runPairedEnd == 0)){
		print SH "\n# Align fastqs to transcriptome for RSEM\n";
		print SH "module load bowtie/1.1.2 \n";
		print SH "export PATH=\$PATH:$NGSbartom/tools/RSEM-1.2.30/\n";
		print SH "date\n";
		print SH "\n# First prepare fastqs\n";
#		print SH "mkdir $outputDirectory/$project/fastq/\n";
		my $fastqstring = $fastqs{$sample};
		$fastqstring =~ s/\,/ /g;
		my $strandstring = "";
#		if ($stranded == 1) { $strandstring = "--strand-specific";}
		if ($stranded == 1) { $strandstring = "--forward-prob 0";}
		#		print SH "gzip -dc $fastqstring | cat > $outputDirectory/$project/fastq/$sample.fastq\n";
		print SH "gunzip $fastqstring\n";
		$fastqstring =~ s/.gz//g;
		print SH "echo $fastqstring\n";
		print SH "cat $fastqstring > $outputDirectory/$project/fastq/$sample.fastq\n";
		print SH "ls $outputDirectory/$project/fastq/$sample.fastq\n";
		print SH "gzip $fastqstring &\n";
		print SH "\n# Then calculate expression for $sample.\n";        
		print SH "rsem-calculate-expression $outputDirectory/$project/fastq/$sample.fastq $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample --no-bam-output -p $numProcessors $strandstring >& $outputDirectory/$project/bam/$sample.rsem.log\n";
		print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n";
		print SH "ls $outputDirectory/$project/fastq/$sample.fastq\n";
		print SH "# rm $outputDirectory/$project/fastq/$sample.fastq\n";		
	    }
	    if ($htseq == 1) {
	    	print SH "module unload mpi\n";
	    	print SH "module load python/anaconda\n";
		print SH "\n# Run htseq-count for $sample\n";
		if ($stranded == 1){
		    print SH "htseq-count -f bam -q -m intersection-nonempty -s reverse -t exon -i gene_id $outputDirectory\/$project\/bam\/$sample.bam $gff{$reference{$sample}} > $outputDirectory\/$project\/bam\/$sample.htseq.counts\n";
		} elsif ($stranded == 0){
		    print SH "htseq-count -f bam -q -m intersection-nonempty -s no -t exon -i gene_id $outputDirectory\/$project\/bam\/$sample.bam $gff{$reference{$sample}} > $outputDirectory\/$project\/bam\/$sample.htseq.counts\n";
		}
	    	print SH "\nmodule unload python/anaconda\n";
	    	print SH "module load gcc/4.8.3\n";
	    }
	    if ($bedtools == 1) {
		print SH "module load bedtools/2.17.0\n";
		print SH "\n# Run bedtools for $sample\n";
	    	print SH "bedtools bamtobed -i $outputDirectory\/$project\/bam\/$sample.bam > $outputDirectory\/$project\/bam\/$sample.bed\n";
		if ($stranded == 1){
		    print SH "bedtools intersect -a $exonbed{$reference{$sample}} -b $outputDirectory\/$project\/bam\/$sample.bed -c -S > $outputDirectory\/$project\/bam\/$sample.bedtools.counts\n";
		} elsif ($stranded == 0) {
		    print SH "bedtools intersect -a $exonbed{$reference{$sample}} -b $outputDirectory\/$project\/bam\/$sample.bed -c > $outputDirectory\/$project\/bam\/$sample.bedtools.counts\n";
		}
		print SH "\n# Clean up bed files.\n";
	    	print SH "rm $outputDirectory\/$project\/bam\/$sample.bed\n";
	    	print SH "\nmodule unload bedtools/2.17.0\n";
	    }
	    print SH "date\n\n";
#	    print SH "samtools flagstat $outputDirectory\/$project\/bam\/$sample.bam > $outputDirectory\/$project\/bam\/$sample.flagstats.txt\n";
	    #	    print SH "date\n";
	    if ($makeTracks == 1){
		print SH "# Create RNA seq Tracks\n";
		print SH "module load R/3.2.2\n";
		if ($stranded == 1){
		    print SH "Rscript $NGSbartom/tools/createRNAseqTracks2.R --assembly=$reference{$sample} --bamDir=$outputDirectory\/$project\/bam\/ --sample=$sample --multiMap=$multiMap\n";
		} else {
		    # The multi mapping argument is used here, because these are mapped with tophat.
		    print SH "Rscript $NGSbartom/tools/createChIPtracks2.R --assembly=$reference{$sample} --bamDir=$outputDirectory\/$project\/bam\/ --sample=$sample --extLen=0 --multiMap=$multiMap\n";
		}
		print SH "date\n\n";
		print SH "mkdir $outputDirectory\/$project\/tracks\n";
		print SH "\n# Make Headers for UCSC genome browser.\n";
		if ($stranded == 1){
		    if ($multiMap == 0) {
			print SH "echo \"track type=bigWig name=$sample.plus.bw description=$sample.plus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=255,0,0 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.plus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.plus.bw.header.txt\n";
			print SH "echo \"track type=bigWig name=$sample.minus.bw description=$sample.minus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.minus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.minus.bw.header.txt\n";
		    } elsif ($multiMap == 1){
			print SH "echo \"track type=bigWig name=$sample.plus.multi.bw description=$sample.plus.multi.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=255,0,0 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.plus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.plus.bw.header.multi.txt\n";
			print SH "echo \"track type=bigWig name=$sample.minus.multi.bw description=$sample.minus.multi.rpm maxHeightPixels=128:60:11 graphtype=bar visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.minus.multi.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.minus.bw.header.multi.txt\n";
		    }
		} else {
		    if ($multiMap == 0){
			print SH "echo \"track type=bigWig name=$sample.bw description=$sample.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.txt\n";
		    }elsif ($multiMap == 1) {
			print SH "echo \"track type=bigWig name=$sample.multi.bw description=$sample.multi.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.multi.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.multi.txt\n";
		    }
		}
		print SH "date\n";
		if ($uploadASHtracks == 1){
		    print SH "# Move tracks into Shilatifard directory structure.\n";
		    print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
		    print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
		    print SH "cp $outputDirectory\/$project\/bam\/$sample*.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
		    if ($uploadBAM == 1){
			print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "# To load, paste the following url into the custom tracks field:\n";
			print SH "# http://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/$scientist/$scientist.$project/$sample.bam\n";
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://$s3path/$scientist.$project/ --region us-west-2\n";
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://$s3path/$scientist.$project/ --region us-west-2\n";
		    }
		    print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
		    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
		    if ($stranded == 1){
			if ($multiMap == 0){
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
		    } else {
			if ($multiMap == 0){
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}			    
		    }			    
		} elsif ($uploadPulmtracks == 1) {
		    if ($uploadBAM == 1){
			print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "# To load, paste the following url into the custom tracks field:\n";
			print SH "# http://s3-us-west-2.amazonaws.com/m-328-data/$scientist/$scientist.$project/$sample.bam\n";
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
		    }
		    print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
		    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
		    if ($stranded == 1){
			if ($multiMap == 0){
			    print SH "aws s3 cp $outputDirectory/bam/$sample.minus.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/bam/$sample.plus.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    print SH "aws s3 cp $outputDirectory/bam/$sample.minus.multi.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/bam/$sample.plus.multi.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}
		    } else {
			if ($multiMap == 0){
			    print SH "aws s3 cp $outputDirectory/bam/$sample.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    print SH "aws s3 cp $outputDirectory/bam/$sample.multi.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}			    
		    }
		    print SH "mv $outputDirectory\/$project\/bam\/$sample.minus.bw $outputDirectory\/$project\/tracks\/\n";
		    print SH "mv $outputDirectory\/$project\/bam\/$sample.plus.bw $outputDirectory\/$project\/tracks\/\n";
		} else {
		    print SH "mv $outputDirectory\/$project\/bam\/$sample.minus.bw $outputDirectory\/$project\/tracks\/\n";
		    print SH "mv $outputDirectory\/$project\/bam\/$sample.plus.bw $outputDirectory\/$project\/tracks\/\n";
		}
	    }
	    close(SH);
	}
    }
    if ($runAlign == 0){
	# Print tips on running the tophat shell scripts.
	print STDERR "To execute all scripts, use the following command:\n";
	print STDERR "find $outputDirectory/*/scripts/ -iname \"*.sh\" -exec msub {} ./ \\\;\n";
    }
}

my %viewpoint;
if ($type eq "4C"){
    open (IN,$fourCdescription);
    while(<IN>) {
	if ($_ !~ /^Sample/){
	    my ($sample,$index,$primer,$assembly,$viewpoint) = split(/\t/,$_);
	    $viewpoint{$sample} = $viewpoint;
	}
    }
}


# If the aligner is bowtie, create the shell scripts to run bowtie on all fastqs, one sample at a time.
if (($buildAlign == 1) && ($aligner eq "bowtie")){
    &datePrint("Creating Bowtie Alignment shell scripts.");
    my @samples;
    # Foreach project (TANGO/MOLNG):
    foreach my $project (keys(%samples)){
	if ($startFromBAM == 0){ $bamDirectory= "$outputDirectory\/$project\/bam";}
	if ($scientists{$project} ne ""){ 
	    $scientist = $scientists{$project};
#	    if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	    $s3path = "ash-tracks/TANGO/$scientist";
	}
	@samples = uniq(split(/\,/,$samples{$project}));
	print STDERR "Samples: @samples\n";
	# Foreach sample within the project:
	foreach my $sample (@samples){
	    # Create a shell script to run bowtie on all fastqs for the sample at the same time.
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_align.sh";
	    &datePrint("Printing $shScript");
	    open (SH,">$shScript");
	    print SH "$header";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    print SH "#MSUB -N $sample\_bowtie\n";
	    print SH "module load bowtie/1.1.2\n";
	    print SH "module load samtools/1.2\n";
	    print SH "module load R/3.2.2\n";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    if ($runTrim == 1){
		print SH "module load java/jdk1.8.0_25\n";
		foreach my $fastq (@fastqs){
		    print SH "\n# Trim poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $fastq.trimmed $trimString\n";
		    print SH "gzip $fastq.trimmed\n";
		    if ($fastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.gz$/){
			if (-f "$1\_fastqc.html"){
			    print SH "\# FastQC file already exists\n";
			} else {
			    print SH "date\n$NGSbartom/tools/FastQC/fastqc $fastq $fastq.trimmed.gz\n";
			}
		    }
		    print SH "mv $fastq.trimmed.gz $fastq\n";
		    print SH "date\n";
		}
	    }
	    my $fastqs = $fastqs{$sample};
	    $fastqs =~ s/,/ /g;
	    print SH "date\n\n";
#	    print STDERR "$sample\tREF:$reference{$sample}\tINDEX:$bowtieIndex{$reference{$sample}}\n";
	    print SH "# Align fastqs with Bowtie\n";
	    if ($multiMap == 0){
		print SH "\ngunzip -c $fastqs | bowtie -p $numProcessors -m 1 -v 2 -S $bowtieIndex{$reference{$sample}} 2> $outputDirectory\/$project\/bam\/$sample.bowtie.log - | samtools view -bS - > $outputDirectory\/$project\/bam\/$sample.bam \n";
	    } elsif ($multiMap ==1) {
		print SH "\ngunzip -c $fastqs | bowtie -p $numProcessors -v 2 -S $bowtieIndex{$reference{$sample}} 2> $outputDirectory\/$project\/bam\/$sample.bowtie.log - | samtools view -bS - > $outputDirectory\/$project\/bam\/$sample.bam \n";
	    }
	    print SH "date\n\n";
	    print SH "# Sort and rearrange bam files\n";
	    print SH "samtools sort $outputDirectory\/$project\/bam\/$sample.bam $outputDirectory\/$project\/bam\/$sample.sorted\n";
	    print SH "date\n";
	    print SH "mv $outputDirectory\/$project\/bam\/$sample.sorted.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n";
	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n\n";
	    if ($type eq "chipseq"){
		if ($makeTracks == 1){
		    print SH "# Make ChIPseq tracks.\n";
		    # The multi mapping argument is not used right now.
		    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
		    print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample --extLen=150\n";
		    print SH "date\n\n";
		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
		    if ($uploadASHtracks == 1){
			print SH "# Move tracks into Shilatifard directory structure\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
			print SH "cp $outputDirectory\/$project\/bam\/$sample*.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			if ($uploadBAM == 1){
			    print SH "cp $bamDirectory\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
		    }
		    if ($uploadPulmtracks == 1){
			if ($uploadBAM == 1){
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/m-328-data/$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}
			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "cp $bamDirectory\/$sample.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
		    }		
		    print SH "mv $bamDirectory\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
		    print SH "\n# Check if bwlist file exists, and if not, create it.\n";
		    print SH "if [ $outputDirectory\/$project\/tracks\/bwlist.txt does not exist ];\nthen\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/bwlist.txt \n";
		    print SH "else\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat >> $outputDirectory\/$project\/tracks\/bwlist.txt \nfi\n";
		    print SH "\n# Make header files for tracks.\n";
		    print SH "echo \"track type=bigWig name=$sample.bw description=$sample.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.txt\n";
		    print SH "date\n\n";
		}
	    }
	    if ($type eq "4C"){
		if ($makeTracks == 1){
		    if (exists($viewpoint{$sample})){
			my ($chr,$start,$stop);
			if ($viewpoint{$sample} =~ /^(chr\w+)\:(\d+)\-(\d+)/){
			    $chr = $1;
			    $start = $2 - 2000;
			    $stop = $3 + 2000;
			    print SH "\n# Remove Viewpoint from Bam files.\n";
			    print SH "samtools view -h $bamDirectory\/$sample.bam | awk '!(\$3 == \"$chr\" && \$4 > $start && \$4 < $stop){print \$0}' | samtools view -Sb - > $bamDirectory\/$sample.noVP.bam\n"; 
			    print SH "\n# Make 4C tracks.\n";
			    # The multi mapping argument is not used right now.
			    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
			    print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample.noVP --extLen=0\n";
			    print SH "mv $bamDirectory\/$sample.noVP.bw $bamDirectory\/$sample.bw\n";
			    print SH "date\n";
			}
		    } else {
			print SH "\n# Make 4C tracks.\n";
			# The multi mapping argument is not used right now.
			# This is because Bowtie doesn't fill in the NH tag in the BAM file.
			print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample --extLen=0\n";
		    }
		    print SH "date\n\n";
		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
		    if ($uploadASHtracks == 1){
			print SH "\n# Move tracks into Shilatifard directory structure\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
			print SH "cp $bamDirectory\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
		    }		
		    print SH "mv $bamDirectory\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
		    print SH "\n# Make header files for tracks.\n";
		    print SH "echo \"track type=bigWig name=$sample.bw description=$sample.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.txt\n";
		    print SH "date\n\n";
		}
	    }
	    close(SH);
	}
    }
    if ($runAlign == 0){
	# Print tips on running the bowtie shell scripts.
	print STDERR "To execute all alignment scripts, use the following command:\n";
	print STDERR "find $outputDirectory/*/scripts/ -iname \"*align.sh\" -exec msub {} ./ \\\;\n";
    }
}


# If the aligner is bowtie, create the shell scripts to run bowtie on all fastqs, one sample at a time.
if (($buildAlign == 1) && ($aligner eq "bwa")){
    &datePrint("Creating BWA Alignment shell scripts.");
    my @samples;
    # Foreach project (TANGO/MOLNG):
    foreach my $project (keys(%samples)){
	if ($startFromBAM == 0){ $bamDirectory= "$outputDirectory\/$project\/bam";}
	if ($scientists{$project} ne ""){ 
	    $scientist = $scientists{$project};
#	    if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	    $s3path = "ash-tracks/TANGO/$scientist";
	}
	@samples = uniq(split(/\,/,$samples{$project}));
	print STDERR "Samples: @samples\n";
	# Foreach sample within the project:
	foreach my $sample (@samples){
	    # Create a shell script to run bwa on all fastqs for the sample at the same time.
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_align.sh";
	    &datePrint("Printing $shScript");
	    open (SH,">$shScript");
	    print SH "$header";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    print SH "#MSUB -N $sample\_bwa\n";
	    print SH "module load bwa/0.7.12\n";
	    print SH "module load samtools/1.2\n";
	    print SH "module load R/3.2.2\n";
	    print SH "module load picard/1.131\n";
	    print SH "module load java/jdk1.8.0_25\n";
	    print SH "\nmkdir $outputDirectory\/$project\/fastq\n";
	    print SH "\nmkdir $outputDirectory\/$project\/fastqc\n";
	    my $PICARD = "/software/picard/1.131/picard-tools-1.131/picard.jar";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    my @newFastqs;
	    if ($runTrim == 1){
		print SH "module load java/jdk1.8.0_25\n";
		foreach my $fastq (@fastqs){
		    print SH "\n# Copy fastq files into outputDirectory.\n";
		    print SH "cp $fastq $outputDirectory/$project/fastq/\n";
		    if ($fastq =~ /\/?([\d\_\-\w\.\.]+)\.fastq\.gz$/){
			$fastq = "$outputDirectory/$project/fastq/$1.fastq.gz";
			push(@newFastqs,$fastq);
		    }
		    print SH "\n# Trim poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $fastq.trimmed $trimString\n";
		    print SH "gzip $fastq.trimmed\n";
		    if ($fastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.gz$/){
			if (-f "$1\_fastqc.html"){
			    print SH "\# FastQC file already exists\n";
			} else {
			    print SH "date\n$NGSbartom/tools/FastQC/fastqc $fastq $fastq.trimmed.gz\n";
			    print SH "mv $fastq*fastqc* $outputDirectory/$project/fastqc/\n";
			}
		    }
		    print SH "mv $fastq.trimmed.gz $fastq\n";
		    print SH "date\n";
		}
		@fastqs = @newFastqs;
	    }
	    my @bams;
	    my $prefix = "";
	    foreach my $fastq (@fastqs){
		print SH "\n# Copy fastq files into outputDirectory.\n";
		print SH "cp $fastq $outputDirectory/$project/fastq/\n";
		if ($fastq =~ /\/?([\d\_\-\w\.]+).fastq\.gz$/){
		    $prefix = $1;
		} else { $prefix = $fastq;}
		$fastq = "$outputDirectory\/$project/fastq/$prefix.fastq.gz";
		push(@newFastqs,$fastq);
	    }
	    @fastqs = @newFastqs;
	    if ($runPairedEnd == 0){
		foreach my $fastq (@fastqs){
		    if ($fastq =~ /\/?([\d\_\-\w\.]+).fastq\.gz$/){
			$prefix = $1;
		    } else { $prefix = $fastq;}
		    print SH "\n# Align Fastq with BWA\n";
		    print SH "bwa mem $bwaIndex{$reference{$sample}} $fastq | samtools view -bS - > $outputDirectory\/$project\/bam\/$prefix.bam\n";
		    print SH "date\n\n";
		    push (@bams,"$outputDirectory\/$project\/bam\/$prefix.bam");
		}
	    } elsif ($runPairedEnd == 1){
		my @read1fastqs;
		my @read2fastqs;
		my @read3fastqs;
		#my @fastqs = split(/\,/,$fastqs{$sample});
		foreach my $fastq (@fastqs){
		    if ($fastq =~ /\_R1\_?\.?/){
			push (@read1fastqs,$fastq);
		    } elsif ($fastq =~ /\_R2\_?\.?/){
			push (@read2fastqs,$fastq);
		    } elsif ($fastq =~ /\_R3\_?\.?/){
			push (@read3fastqs,$fastq);
		    } else {
			print STDERR "ERR: Could not find Read Number!\n";
		    }
		}
		@read1fastqs = sort(@read1fastqs);
		@read2fastqs = sort(@read2fastqs);
		@read3fastqs = sort(@read3fastqs);
		my $read1fastqs = "@read1fastqs";
		my $read2fastqs = "@read2fastqs";
		my $read3fastqs = "@read3fastqs";
		if (length($read3fastqs)>length($read2fastqs)){
		    @read2fastqs = @read3fastqs;
		}
		$read1fastqs =~ s/\ /,/g;
		$read2fastqs =~ s/\ /,/g;
		for (my $i=0;$i <= $#read1fastqs;$i++){
		    my $prefix = "";
		    if ($read1fastqs[$i] =~ /\/?([\d\_\-\w\.]+)\_R1\.fastq\.gz$/){
			$prefix = $1;
		    } else { $prefix = $read1fastqs[$i];}
		    print SH "\n# Align Fastq with BWA\n";
		    print SH "bwa mem $bwaIndex{$reference{$sample}} $read1fastqs[$i] $read2fastqs[$i] | samtools view -bS - > $outputDirectory\/$project\/bam/$prefix.bam\n";
		    print SH "date\n\n";
		    push (@bams,"$outputDirectory\/$project\/bam\/$prefix.bam");
		}
	    }
	    my $bamString = "I=@bams";
	    $bamString =~ s/\s/ I=/g;
	    print SH "# Merge bam files from the same sample $sample\n";
	    print SH "java -Xmx2g -jar $PICARD MergeSamFiles $bamString O=$outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n\n";
	    print SH "# Mark duplicate reads in the bam file for sample $sample\n";
	    print SH "java -Xmx2g -jar $PICARD MarkDuplicates I=$outputDirectory\/$project\/bam\/$sample.bam O=$outputDirectory\/$project\/bam\/$sample.mdup.bam M=$outputDirectory\/$project\/bam\/$sample.mdup_metrics.txt\n";
 	    print SH "date\n\n";
	    print SH "# Replace unmarked bam with marked bam.\n";
 	    print SH "mv $outputDirectory\/$project\/bam\/$sample.mdup.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
 	    print SH "date\n";
	    print SH "# Index bam file and gather flagstats.\n";
 	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "samtools flagstat $outputDirectory\/$project\/bam\/$sample.bam > $outputDirectory\/$project\/bam\/$sample.bam.flagstats.txt\n";
 	    print SH "date\n\n";
 	    if ($type eq "chipseq"){
 		if ($makeTracks == 1){
 		    print SH "# Make ChIPseq tracks.\n";
 		    # The multi mapping argument is not used right now.
 		    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
 		    print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample --extLen=150\n";
 		    print SH "date\n\n";
 		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
 		    if ($uploadASHtracks == 1){
 			print SH "# Move tracks into Shilatifard directory structure\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
 			print SH "cp $outputDirectory\/$project\/bam\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			if ($uploadBAM == 1){
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
 		    }
		    if ($uploadPulmtracks == 1){
			if ($uploadBAM == 1){
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/m-328-data$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp $outputDirectory\/$project\/bam\/$sample.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
 		    }		
 		    print SH "mv $outputDirectory\/$project\/bam\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
 		    print SH "\n# Check if bwlist file exists, and if not, create it.\n";
 		    print SH "if [ $outputDirectory\/$project\/tracks\/bwlist.txt does not exist ];\nthen\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/bwlist.txt \n";
 		    print SH "else\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat >> $outputDirectory\/$project\/tracks\/bwlist.txt \nfi\n";
 		    print SH "\n# Make header files for tracks.\n";
 		    print SH "echo \"track type=bigWig name=$sample.bw description=$sample.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.txt\n";
 		    print SH "date\n\n";
 		}
 	    }
 	    if ($type eq "4C"){
 		if ($makeTracks == 1){
 		    if (exists($viewpoint{$sample})){
 			my ($chr,$start,$stop);
 			if ($viewpoint{$sample} =~ /^(chr\w+)\:(\d+)\-(\d+)/){
 			    $chr = $1;
 			    $start = $2 - 2000;
 			    $stop = $3 + 2000;
 			    print SH "\n# Remove Viewpoint from Bam files.\n";
 			    print SH "samtools view -h $bamDirectory\/$sample.bam | awk '!(\$3 == \"$chr\" && \$4 > $start && \$4 < $stop){print \$0}' | samtools view -Sb - > $bamDirectory\/$sample.noVP.bam\n"; 
 			    print SH "\n# Make 4C tracks.\n";
 			    # The multi mapping argument is not used right now.
 			    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
 			    print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample.noVP --extLen=0\n";
 			    print SH "mv $bamDirectory\/$sample.noVP.bw $bamDirectory\/$sample.bw\n";
 			    print SH "date\n";
 			}
 		    } else {
 			print SH "\n# Make 4C tracks.\n";
 			# The multi mapping argument is not used right now.
 			# This is because Bowtie doesn't fill in the NH tag in the BAM file.
 			print SH "Rscript $NGSbartom/tools/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample --extLen=0\n";
 		    }
 		    print SH "date\n\n";
 		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
 		    if ($uploadASHtracks == 1){
 			print SH "\n# Move tracks into Shilatifard directory structure\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
 			print SH "cp $outputDirectory\/$project\/bam\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
 		    }		
 		    print SH "mv $outputDirectory\/$project\/bam\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
 		    print SH "\n# Make header files for tracks.\n";
 		    print SH "echo \"track type=bigWig name=$sample.bw description=$sample.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.bw.header.txt\n";
 		    print SH "date\n\n";
		}
	    }
	    close(SH);
	}
	     }
     if ($runAlign == 0){
 	# Print tips on running the bwa shell scripts.
 	print STDERR "To execute all alignment scripts, use the following command:\n";
 	print STDERR "find $outputDirectory/*/scripts/ -iname \"*align.sh\" -exec msub {} ./ \\\;\n";
     }
}


if (($buildAlign ==1) && ($runAlign ==1)){
    # Submit the alignment jobs, saving the job id.
    # This will finish each project (TANGO) before starting the next.  Is this what we want?
    foreach my $project (keys(%samples)){
	if ($scientists{$project} ne ""){ 
	    $scientist = $scientists{$project};
#	    if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	    $s3path = "ash-tracks/TANGO/$scientist";
	}
	&datePrint ("Starting alignment scripts.");
	my $result = `find $outputDirectory/*/scripts/ -iname \"*align.sh\" -exec msub {} ./ \\\;`;
	$result =~ s/\s+/\:/g;
	$result =~ s/^://;
	$result =~ s/:$//;
	&datePrint( "Need to wait for the following jobs to finish:\n$result");
	open (SH, ">$outputDirectory/$project/scripts/AlignmentDependentScript.sh");
	print SH $header;
	print SH "#MSUB -W depend=afterok:$result\n";
	print SH "#MSUB -N PostAlignmentAnalysis\n";
	print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	print SH "\necho \"Alignment jobs $result have finished.\"\n";
	if ($aligner eq "tophat"){
	    print SH "module load R/3.2.2\n";
	    print SH "\n# Make Tophat report summarizing alignment.\n";
	    print SH "Rscript $NGSbartom/tools/createTophatReport.R --topHatDir=$outputDirectory\/$project\/Tophat_aln --nClus=$numProcessors\n";
	}elsif ($aligner eq "bowtie"){
	    print SH "find $bamDirectory\/*bowtie.log -type f -print -exec cat {} \; >> $outputDirectory\/$project\/alignlog.txt\n";
	}
	close SH;
	&datePrint("Creating dependent job that will only run after alignments finish.");
	my $result2 = `msub $outputDirectory/$project/scripts/AlignmentDependentScript.sh`;
	$result2 =~ s/\s+//g;
	my $jobfinished = "no";
	# Wait until the job is Complete.
	&datePrint("Waiting for job $result2 to finish. (each . = 30 seconds)");
	# Check qstat every 30 seconds, adding a "." every time you check.
	until ($jobfinished eq "Completed"){
	    $jobfinished = `checkjob $result2 | grep ^State:`;
	    sleep(30);
	    print STDERR ".";
	    if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
	    print STDERR "$jobfinished";
	}
	print STDERR "\n";
	&datePrint("Job $result2 done.  Continuing.");
    }
}

if ($buildGenotyping ==1) {
    if ($type eq "RNA"){
	foreach my $project (keys(%samples)){
	    my @samples = uniq(split(/\,/,$samples{$project}));
	    if (-d "$outputDirectory\/$project\/genotype/") {
		print STDERR "Genotype Directory found.\n";
	    }else {system("mkdir $outputDirectory\/$project\/genotype");print STDERR "Genotype directory created.\n";}
	    # Foreach sample within the project:		    
	    foreach my $sample (@samples){
		open (SH, ">$outputDirectory/$project/scripts/$sample\_genotype.sh");
		print SH $header;
		print SH "#MSUB -N $sample\_genotype\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		print SH "module load samtools/1.2\n";
		print SH "\n\n";
		print SH "# Sort BAM file.\n";
		print SH "$NGSbartom/tools/samtools-0.1.19/samtools sort $outputDirectory\/$project\/bam\/$sample.bam $outputDirectory\/$project\/bam\/$sample.sorted\n";
		print SH "date\n\n";
		print SH "# Mark Duplicates with Picard.\n";
		print SH "java -jar $NGSbartom/tools/picard.jar MarkDuplicates I=$outputDirectory\/$project\/bam\/$sample.sorted.bam O=$outputDirectory\/$project\/bam\/$sample.mdup.bam M=$outputDirectory\/$project\/bam\/$sample.mdup.metrics.txt\n";
		print SH "date\n\n";
#		print SH "# Sort mdup BAM file with Picard.\n";
#		print SH "java -jar $NGSbartom/tools/picard.jar SortSam I=$outputDirectory\/$project\/bam\/$sample.mdup.bam O=$outputDirectory\/$project\/bam\/$sample.mdup.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true\n";
		#		print SH "date\n\n";
		print SH "# Reorder mdup BAM file with Picard.\n";
		print SH "java -jar $NGSbartom/tools/picard.jar ReorderSam I=$outputDirectory\/$project\/bam\/$sample.mdup.bam O=$outputDirectory\/$project\/bam\/$sample.mdup.reordered.bam R=$gatkRef{$reference{$sample}} CREATE_INDEX=true\n";
		print SH "date\n\n";
		print SH "# Split Reads at splicing events (runs of Ns in CIGAR string)\n";
		print SH "java -jar $NGSbartom/tools/GATK_v3.6/GenomeAnalysisTK.jar -T SplitNCigarReads -R $gatkRef{$reference{$sample}} -I $outputDirectory\/$project\/bam\/$sample.mdup.reordered.bam -o $outputDirectory\/$project\/bam\/$sample.split.bam -U ALLOW_N_CIGAR_READS -fixNDN\n";
		print SH "date\n\n";
		print SH "# Find Target regions for Realignment.\n";
		print SH "java -jar $NGSbartom/tools/GATK_v3.6/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $gatkRef{$reference{$sample}} -I $outputDirectory\/$project\/bam\/$sample.split.bam -o $outputDirectory\/$project\/bam\/$sample.split.intervals.list --known $knownIndelsites{$reference{$sample}}\n";
		print SH "date\n\n";
		print SH "# Realign indels in target regions.\n";
		print SH "java -jar $NGSbartom/tools/GATK_v3.6/GenomeAnalysisTK.jar -T IndelRealigner -R $gatkRef{$reference{$sample}} -I $outputDirectory\/$project\/bam\/$sample.split.bam -targetIntervals $outputDirectory\/$project\/bam\/$sample.split.intervals.list -known $knownIndelsites{$reference{$sample}} -o $outputDirectory\/$project\/bam\/$sample.split.real.bam\n";
		print SH "date\n\n";
		print SH "# Generating Base Recalibration Table.\n";
		print SH "java -jar $NGSbartom/tools/GATK_v3.6/GenomeAnalysisTK.jar -T BaseRecalibrator -R $gatkRef{$reference{$sample}} -I $outputDirectory\/$project\/bam\/$sample.split.real.bam -o $outputDirectory\/$project\/genotype\/$sample.split.real.recal.table -knownSites $knownSNPsites{$reference{$sample}}\n";
		print SH "date\n\n";
		print SH "# Calling SNPs and Indels with HaplotypeCaller.\n";
		print SH "java -jar $NGSbartom/tools/GATK_v3.6/GenomeAnalysisTK.jar -T HaplotypeCaller -R $gatkRef{$reference{$sample}} -I $outputDirectory\/$project\/bam\/$sample.split.real.bam -o $outputDirectory\/$project\/genotype\/$sample.raw.snps.indels.vcf --dbsnp $knownSNPsites{$reference{$sample}}\n";
		print SH "date\n\n";
		close SH;
	    }
	    if ($runGenotyping == 1){
		&datePrint("Submitting jobs for genotype analysis.");
		my $result = `find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec msub {} ./ \\\;`;
	    } else {
		print STDERR "Not submitting genotype scripts.  To run them use: \`find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec msub \{\} ./ \\\;\`\n";
	    }
	}
    } else { print "ERR: Genotyping not yet implemented for whole genome or exome sequencing.\n";}
}


if ($buildEdgeR ==1) {
    foreach my $project (keys(%samples)){
	my @countsfiles;
	open (SH, ">$outputDirectory/$project/scripts/downstreamRNAanalysis.sh");
	print SH $header;
	print SH "#MSUB -N downstreamRNAanalysis\n";
	print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	print SH "module load R/3.2.2\n";
	if ($type eq "RNA"){
	    if ($startFromBAM == 1){
		# If the bam directory was specified in the input (no base space directory or fastq directory)
		if ($htseq == 1 || $bedtools == 1) {
		    my @samples = uniq(split(/\,/,$samples{$project}));
		    # Foreach sample within the project:		    
		    if ($htseq == 1) {
		    	print SH "module unload mpi\n";
	    		print SH "module load python/anaconda\n";
	    		my $htseq_sample_count = 0;
			print SH "\n# Launch up to $numProcessors htseq jobs\n";
	    		foreach my $sample (@samples){
			    if ($stranded == 1){
				print SH "htseq-count -f bam -q -m intersection-nonempty -s reverse -t exon -i gene_id $bamDirectory\/$sample.bam $gff{$reference{$sample}} > $bamDirectory\/$sample.htseq.counts &\n";
			    } elsif ($stranded == 0){
				print SH "htseq-count -f bam -q -m intersection-nonempty -s no -t exon -i gene_id $bamDirectory\/$sample.bam $gff{$reference{$sample}} > $bamDirectory\/$sample.htseq.counts &\n";
			    }
			    $htseq_sample_count++;
			    if ($htseq_sample_count % $numProcessors == 0) {
				print SH "\n# Wait for htseq jobs to finish, then start more.\n";
				print SH "wait\n";
			    }
	    		}
	    		if ($htseq_sample_count % $numProcessors > 0) {
			    print SH "wait\n";
			}    			
	    		print SH "module unload python/anaconda\n";
	    		print SH "module load gcc/4.8.3\n";
		    }
		    if ($bedtools == 1) {
	    		print SH "module load bedtools/2.17.0\n";
	    		my $bedtools_sample_count = 0;
			print SH "\n# Launch up to $numProcessors bamtobed jobs\n";
	    		foreach my $sample (@samples){
			    print SH "bedtools bamtobed -i $bamDirectory\/$sample.bam > $bamDirectory\/$sample.bed &\n";
			    $bedtools_sample_count++;
			    if ($bedtools_sample_count % $numProcessors == 0) {
				print SH "\n# Wait for bedtools jobs to finish.\n";
				print SH "wait\n";
			    }
	    		}
	    		if ($bedtools_sample_count % $numProcessors > 0) {
			    print SH "\n# Wait for bedtools jobs to finish, then start more.\n";
			    print SH "wait\n";
	    		}
	    		my $bedtools_sample_count2 = 0;
			print SH "\n# Launch up to $numProcessors bedtools intersect jobs\n";
	    		foreach my $sample (@samples){
			    if ($stranded == 1){
				print SH "bedtools intersect -a $exonbed{$reference{$sample}} -b $bamDirectory\/$sample.bed -c -S > $bamDirectory\/$sample.bedtools.counts &\n";
			    } elsif ($stranded == 0) {
				print SH "bedtools intersect -a $exonbed{$reference{$sample}} -b $bamDirectory\/$sample.bed -c > $bamDirectory\/$sample.bedtools.counts &\n";
			    }
			    $bedtools_sample_count2++;
			    if ($bedtools_sample_count2 % $numProcessors == 0) {
				print SH "wait\n\n";
			    }
	    		}
	    		if ($bedtools_sample_count2 % $numProcessors > 0) {
			    print SH "wait\n\n";
	    		}
			print SH "\n# Clean up bed files.\n";
	    		foreach my $sample (@samples){
			    print SH "rm $bamDirectory\/$sample.bed\n";
	    		}
	    		print SH "\nmodule unload bedtools/2.17.0\n";
		    }
		}
		print SH "module load R/3.2.2\n";	       				
		if ($makeTracks == 1){
		    my @samples = uniq(split(/\,/,$samples{$project}));
		    # Foreach sample within the project:
		    foreach my $sample (@samples){
		    # if the user wants tracks, put that in the downstream analysis script.
			print SH "# Create RNA seq Tracks\n";
			print SH "Rscript $NGSbartom/tools/createRNAseqTracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory\/ --sample=$sample\n";
			print SH "date\n\n";
			print SH "mkdir $outputDirectory\/$project\/tracks\n";
			print SH "\n# Make Headers for UCSC genome browser.\n";
			print SH "echo \"track type=bigWig name=$sample.plus.bw description=$sample.plus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=255,0,0 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.plus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.plus.bw.header.txt\n";
			print SH "echo \"track type=bigWig name=$sample.minus.bw description=$sample.minus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.minus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.minus.bw.header.txt\n";
			print SH "date\n";
			if ($uploadASHtracks == 1){
			    print SH "# Move tracks into Shilatifard directory structure.\n";
			    print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
			    print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
			    print SH "cp $bamDirectory\/$sample.minus.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			    print SH "cp $bamDirectory\/$sample.plus.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			    print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";			
			} else {
			    print SH "mv $bamDirectory\/$sample.minus.bw $outputDirectory\/$project\/tracks\/\n";
			    print SH "mv $bamDirectory\/$sample.plus.bw $outputDirectory\/$project\/tracks\/\n";
			}
		    }
		}
	    } else { $bamDirectory = "$outputDirectory\/$project\/bam\/";}
	    print SH "\n# Make Analysis directory for all analysis files.\n";
	    print SH "mkdir $outputDirectory\/$project\/analysis\n";
	    my @methods;
	    if ($granges == 1){
		push (@methods,"granges");
		print SH "\n# Make granges counts table for downstream analysis.\n";
		print SH "Rscript $NGSbartom/tools/createRNAcounts2.R --assembly=$reference{$project} --bamDir=$bamDirectory --numCores=$numProcessors --txdbfile=$txdbfile{$reference{$project}}\n";
		print SH "\n# Rename counts.txt files to specify counting method\n";
		print SH "mv $bamDirectory/counts.txt $bamDirectory/granges.all.counts.txt\n";
		print SH "mv $bamDirectory/counts.rda $bamDirectory/granges.all.counts.rda\n";
		print SH "mv $bamDirectory/gnModel.txt $bamDirectory/granges.all.gnModel.txt\n";
		print SH "mv $bamDirectory/gnModel.rda $bamDirectory/granges.all.gnModel.rda\n";
		print SH "date\n";
		push (@countsfiles,"$bamDirectory/granges.all.counts.txt");
	    }
	    if ($htseq == 1) {
		push (@methods,"htseq");
		print SH "\n# Make HTseq counts table.\n";
	    	print SH "perl $NGSbartom/tools/makeHTseqCountsTable.pl $bamDirectory\/ $gff{$reference{$project}} $bamDirectory\/\n";
		push (@countsfiles,"$bamDirectory/htseq.all.counts.txt");
	    }
	    if ($rsem == 1) {
		push (@methods,"rsem");
		print SH "\n# Make RSEM counts table.\n";
	    	print SH "perl $NGSbartom/tools/makeRSEMcountsTable.pl $bamDirectory\/ $gff{$reference{$project}} $bamDirectory\/ isoforms\n";
		push (@countsfiles,"$bamDirectory/rsem.all.counts.txt");
	    }
	    if ($bedtools == 1) {
		push (@methods,"bedtools");
		print SH "\n# Make Bedtools counts table.\n";
	    	print SH "perl $NGSbartom/tools/makeBEDtoolsCountsTable.pl $bamDirectory\/ $bamDirectory\/\n";
		push (@countsfiles,"$bamDirectory/htseq.all.counts.txt");
	    }
	    if ($comparisons ne ""){
		foreach my $method (@methods){
		    print SH "\n# Run EdgeR, using comparisons file, without MDS plot (which sometimes crashes), $method.\n";
		    print SH "Rscript $NGSbartom/tools/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --comparisonFile=$comparisons --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis --runMDS=0\n";
		    print SH "\n# Run EdgeR, creating MDS plot, but not running comparisons, $method.\n";
		    print SH "Rscript $NGSbartom/tools/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis --runMDS=1\n";
		    print SH "date\n";
		    if ($method eq "granges"){
			print SH "\n# Create Correlation plot for all samples\n";
			print SH "# This plot is still in development, so don't over-interpret, $method.\n";
			print SH "Rscript $NGSbartom/tools/makeCorrelationPlotAllSamples.R --countFile=$bamDirectory\/$method.all.counts.rda --outputDirectory=$outputDirectory\/$project\/analysis\n";
		    }
		}
		print SH "date\n";
		open(CMP,$comparisons);
		while(<CMP>){
		    chomp $_;
		    my ($comp,@groups) = split(/\,/,$_);
		    if ("@groups" =~ /^[1\-0\s]+/){
			foreach my $method (@methods){
			    print SH "\n# Create MA plot for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/makeMAplot.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --labelTop=1\n";
			    print SH "date\n";
			    print SH "\n# Create Big Heatmaps for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/makeBigHeatmap.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --countFile=$outputDirectory\/$project\/analysis\/$method.normCounts.txt\n";
			    print SH "date\n";
			    print SH "\n# Run GO analysis for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/runGOforDEG.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --assembly=$reference{$project}\n";
			    print SH "date\n";
			}
		    }
		}
		close(CMP);
	    } else {
		foreach my $method (@methods){
		    print SH "\n# Create MDS plot for samples, with count method $method.\n";
		    print SH "Rscript $NGSbartom/tools/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis\n";
		}		
	    }
	    if ($runEdgeR == 1){
		&datePrint("Submitting job for downstream RNA-seq analysis.");
		$cmd = "msub $outputDirectory/$project/scripts/downstreamRNAanalysis.sh";
		system($cmd);
	    }
	}
	close SH;
    }
}

if (($buildPeakCaller ==1) && ($type eq "chipseq")){
    foreach my $project (keys(%samples)){
	&datePrint("Initiating peak calling for project $project");
	my $ref = $reference{$project};
	$ref =~ s/\d+//g;
	if ($ref eq "sacCer"){ $ref = "ce";}
	if ($ref eq "hg"){ $ref = "hs";}
	if ($startFromBAM == 0){ $bamDirectory = "$outputDirectory\/$project\/bam\/";}
	if ($chipDescription ne ""){
	    $cmd = "mkdir $outputDirectory\/$project\/peaks";
	    system($cmd);
	    my $broadPeakCheck = `grep -i \"broad\" $chipDescription`;
	    if ($broadPeakCheck ne ""){
		&datePrint("Found some broad peaks to call. Creating SICER script.");
		open (BSH, ">$outputDirectory/$project/scripts/runSICER\_callPeaks.sh");
		print BSH $header;
		print BSH "#MSUB -N callSicerPeaks\n";
		print BSH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		print BSH "\nmodule unload R\n";
		print BSH "module unload mpi\n";
		print BSH "module load python/anaconda\n";
		print BSH "module load bedtools/2.17.0\n";
		print BSH "module load samtools/1.2\n";
		print BSH "export PATH=$NGSbartom/tools/SICER_V1.1/SICER/:\$PATH\n";
	    }
	    open(CHIP,$chipDescription);
	    while(<CHIP>){
		if (($_ !~ /^IP/)&& ($_ !~ /^ChIP/)){
		    chomp $_;
		    my ($ip,$input,$peakType,$sample,$peakset) = split(/\,/,$_);
		    $peakType = lc $peakType;
		    &datePrint("Setting up scripts for $peakType peaks for IP file $ip");
		    if ($peakType eq "narrow"){
			open (SH, ">$outputDirectory/$project/scripts/run\_$ip\_callPeaks.sh");
			print SH $header;
			print SH "#MSUB -N $ip\_NarrowPeaks\n";
			print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
			#print SH "module load R/3.2.2\n";
			print SH "export PATH=\$PATH:$NGSbartom/tools/MACS-1.4.2/bin\n";
			print SH "export PYTHONPATH=$NGSbartom/tools/MACS-1.4.2/lib/python2.6/site-packages:\$PYTHONPATH\n";
			print SH "\n# Call $peakType peaks for ip file $ip, input $input\n";
			print SH "macs14 -t $bamDirectory\/$ip.bam -c $bamDirectory\/$input.bam -f BAM -g $ref -n $outputDirectory\/$project\/peaks\/$ip.macsPeaks >& $outputDirectory\/$project\/peaks\/$ip.macs14.log\n";
			print SH "date\n";
			if ($makeTracks){
			    print SH "\n# Convert bed files to bb files for upload to genome browser.\n";
			    print SH "# First convert float scores to integers, and get rid of any scores over the maximum (1000)\n";
			    print SH "perl -nle \' print \"\$1\\t\$2\\t\$3\\t\$4\\t\$5\" if \/^\(\\w+\)\\t\(\\d+\)\\t\(\\d+\)\\t\(\[\\w\\_\\-\\.\]+\)\\t\(\\d+\)\/ && \(\$5 <= 1000)\' $outputDirectory\/$project\/peaks\/$ip.macsPeaks\_peaks.bed > $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.int.bed\n";
			    print SH "# Then add back in the high-scoring peaks, now with a score of 1000.\n";
			    print SH "perl -nle \' print \"\$1\\t\$2\\t\$3\\t\$4\\t\$5\" if \/^\(\\w+\)\\t\(\\d+\)\\t\(\\d+\)\\t\(\[\\w\\_\\-\\.\]+\)\\t\(\\d+\)\/ \' $outputDirectory\/$project\/peaks\/$ip.macsPeaks\_peaks.bed | awk \'\$5>1000 \{ print \$1,\"\\t\",\$2,\"\\t\",\$3,\"\\t\",\$4,\"\\t1000\"}\' \>\> $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.int.bed\n";
			    print SH "# Sort the bed file to get everything back in chromosomal order.\n";
			    print SH "sort -k1,1 -k2,2n $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.int.bed > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.capped.bed\n";
			    print SH "rm $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.int.bed\n";
			    print SH "perl -pe \"s/ //g\" $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.bed > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed\n";
#			    print SH "ln -s $outputDirectory\/$project\/peaks\/$ip.macsPeaks_peaks.bed $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed\n";
			    print SH "# Now do the conversion with bedToBigBed\n";
			    print SH "$NGSbartom/tools/bedToBigBed $outputDirectory\/$project\/peaks\/$ip.macsPeaks.capped.bed $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb\n";
			    print SH "date\n";
			    if ($uploadASHtracks == 1){
				print SH "\n# Copy bigBed to Amazon S3, for UCSC genome browser to access.\n";
				print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
				print SH "aws s3 cp $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb s3://$s3path/$scientist.$project/ --region us-west-2\n";
				print SH "date\n";
			    }
			    print SH "\n# Create track description for bigBed file.\n";
			    print SH "echo \"track type=bigBed name=$ip.macsPeaks description=\\\"MACS peaks in $ip relative to $input\\\" graphtype=bar maxHeightPixels=128:60:11 visibility=dense color=0,0,0 itemRGB=on useScore=1 autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$ip.macsPeaks.bb\" | cat > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb.header.txt\n";
			    print SH "date\n";
			}
			print SH "\n# Annotate peaks with nearby genes.\n";
			print SH "module load R/3.2.2\n";
			print SH "Rscript $NGSbartom/tools/addGenesToBed.R --peakFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed --outputDirectory=$outputDirectory\/$project\/peaks --assembly=$reference{$project} --txdbfile=$txdbfile{$reference{$project}}\n";
			print SH "\n# Extend peaks from the summit, adding $upstream bp upstream and $downstream bp downstream.\n";
			print SH "module load bedtools/2.17.0\n";
			print SH "bedtools slop -i $outputDirectory\/$project\/peaks\/$ip.macsPeaks\_summits.bed -g $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes -l $upstream -r $downstream > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.bed\n";
			print SH "\n# Filter out peaks with an maximum input rpm over 1.\n";
			print SH "Rscript $NGSbartom/tools/filterOutHighInputPeaks.R  --inputfile=$outputDirectory\/$project\/tracks\/$input.bw --bedfile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.bed --maxInput=1\n";
			print SH "\n# Take the top peaks (at most 5000) and continue with them.\n";
			print SH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.bed | head -n 5000 > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.top5k.bed\n";
			print SH "\n# Make a heatmap and meta plot for expanded peaks.\n";
			print SH "Rscript $NGSbartom/tools/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.top5k.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";

			print SH "\n# Find Top 100 TSS-proximal peaks, find coordinates of regions around associated TSS's  and make a heatmap and meta plot.\n";
			print SH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.macsPeaks.anno.txt | awk \'\$10 < $distToTSS {print \$5,\$10,\$11}\' | awk \'\{print \$3\}\' | sort | uniq | head -n 100 > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.txt\n";
			print SH "Rscript $NGSbartom/tools/fromGeneListToTSSbed.R --txdbfile=$txdbfile{$reference{$project}} --assembly=$reference{$project} --geneList=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.txt --up=$upstream --down=$downstream --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			print SH "Rscript $NGSbartom/tools/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.$upstream.$downstream.tss.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			close SH;
		       
		    } elsif ($peakType eq "broad"){
			print BSH "module load bedtools/2.17.0\n";
			print BSH "\n# Convert bam files to bed files, as needed.\n";
			print BSH "if \[ ! -s \"$bamDirectory\/$ip.bed\" ]; then\n";
			print BSH "\tbedtools bamtobed -i $bamDirectory\/$ip.bam > $bamDirectory\/$ip.bed\n";
			print BSH "\tdate\n";
			print BSH "fi\n";
			print BSH "if \[ ! -s \"$bamDirectory\/$input.bed\" ]; then\n";
			print BSH "\tbedtools bamtobed -i $bamDirectory\/$input.bam > $bamDirectory\/$input.bed\n";
			print BSH "\tdate\n";
			print BSH "fi\n";
			print BSH "module unload bedtools\n";
			print BSH "\n# Call $peakType peaks for ip file $ip, input $input\n";
			print BSH "mkdir $outputDirectory\/$project\/peaks\/$ip.sicer\n";
			print BSH "SICER.sh $bamDirectory $ip.bed $input.bed $outputDirectory\/$project\/peaks\/$ip.sicer $reference{$project} 1 200 150 0.8 600 1e-8 >& $outputDirectory\/$project\/peaks\/$ip.sicer.log\n";
			# $ sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["FDR"]			
			print BSH "date\n";
			print BSH "ln -s $outputDirectory\/$project\/peaks\/$ip.sicer\/$ip-W200-G600-FDR1e-8-island.bed $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed\n";
			if ($makeTracks){
			    print BSH "\n# NB: Changing peak scores for bigBed purposes! Use original file for original scores!\n";
			    print BSH "\n# Remove peaks with scores over 1000.\n";
			    print BSH "awk \'\$4 <= 1000 {printf \"\%s\\t\%d\\t\%d\\t\%s\_\%d\\t\%d\\n\", \$1,\$2,\$3,\"$ip\",NR,\$4}\' $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed\n";
			    print BSH "\n# Then add back in the high-scoring peaks, now with a score of 1000.\n";
			    print BSH "awk \'\$4 > 1000 {printf \"\%s\\t\%d\\t\%d\\t\%s\_\%d\\t\%d\\n\", \$1,\$2,\$3,\"$ip\",NR,1000}\' $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed >> $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed\n";
			    print BSH "\n# Sort the bed file to get everything back in chromosomal order.\n";
			    print BSH "sort -k1,1 -k2,2n $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.sorted.bed\n";
			    print BSH "mv $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.sorted.bed $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed\n";
			    print BSH "\n# Now do the conversion with bedToBigBed\n";
			    print BSH "$NGSbartom/tools/bedToBigBed $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bb\n";
			    print BSH "date\n";
			    if ($uploadASHtracks == 1){
				print BSH "\n# Copy bigBed to Amazon S3, for UCSC genome browser to access.\n";
				print BSH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
				print BSH "aws s3 cp $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bb s3://$s3path/$scientist.$project/ --region us-west-2\n";
				print BSH "date\n";
			    }
			    print BSH "\n# Create track description for bigBed file.\n";
			    print BSH "echo \"track type=bigBed name=$ip.sicerPeaks description=\\\"SICER peaks in $ip relative to $input\\\" graphtype=bar maxHeightPixels=128:60:11 visibility=dense color=0,0,0 itemRGB=on useScore=1 autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$ip.sicerPeaks.bb\" | cat > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bb.header.txt\n";
			    print BSH "date\n";
			}
			print BSH "\n# Add a peak name to all peaks.\n";
			print BSH "awk \'{printf \"\%s\\t\%d\\t\%d\\t\%s\_\%d\\t\%d\\n\", \$1,\$2,\$3,\"$ip\",NR,\$4}\' $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.named.bed\n";
			print BSH "mv $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.named.bed $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed\n";
			print BSH "date\n";
			print BSH "\n# Annotate peaks with nearby genes.\n";
			print BSH "\nmodule unload python\nmodule load mpi/openmpi-1.6.3-gcc-4.6.3\nmodule load R/3.2.2\n";
			print BSH "Rscript $NGSbartom/tools/addGenesToBed.R --peakFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed --outputDirectory=$outputDirectory\/$project\/peaks --assembly=$reference{$project} --txdbfile=$txdbfile{$reference{$project}}\n";
			print BSH "date\n";
			print BSH "\n# Find summits of peaks.\n";
			print BSH "Rscript $NGSbartom/tools/fromBedPlusBWtoSummit.R --bedfile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			print BSH "module load bedtools/2.17.0\n";
			print BSH "bedtools slop -i $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.summits.bed -g $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes -l $upstream -r $downstream > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.bed\n";
			print BSH "\n# Filter out peaks with an maximum input rpm over 1.\n";
			print BSH "Rscript $NGSbartom/tools/filterOutHighInputPeaks.R  --inputfile=$outputDirectory\/$project\/tracks\/$input.bw --bedfile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.bed --maxInput=1\n";
			print BSH "\n# Make a heatmap and meta plot for top 5000 broad peaks.\n";
			print BSH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.bed | head -n 5000 > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.top5k.bed\n";
			print BSH "Rscript $NGSbartom/tools/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.top5k.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			print BSH "\n# Find Top 100 TSS-proximal peaks, find coordinates of regions around associated TSS's  and make a heatmap and meta plot.\n";
			print BSH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.anno.txt | awk \'\$10 < $distToTSS {print \$5,\$10,\$11}\' | awk \'\{print \$3\}\' | sort | uniq | head -n 100 > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.txt\n";
			print BSH "Rscript $NGSbartom/tools/fromGeneListToTSSbed.R --txdbfile=$txdbfile{$reference{$project}} --assembly=$reference{$project} --geneList=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.txt --up=$upstream --down=$downstream --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			print BSH "Rscript $NGSbartom/tools/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.$upstream.$downstream.tss.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			print BSH "\nmodule unload R\nmodule unload mpi\nmodule load python/anaconda\n";
		    }		
		}
	    }
	    close(CHIP);
	} else {
	    die "ERR:  ChipDescription file must be provided which specifies pairs of Input and IP samples.\n";
	}
    }
    if ($runPeakCaller ==1){
	# Submit the peak calling jobs, saving the job id.
	## This will finish each project (TANGO) before starting the next.  Is this what we want?
    	foreach my $project (keys(%samples)){
    	    &datePrint ("Starting peak calling scripts.");
    	    my $result = `find $outputDirectory/*/scripts/ -iname \"*callPeaks.sh\" -exec msub {} ./ \\\;`;
	    $result =~ s/\s+/\:/g;
	    $result =~ s/^://;
	    $result =~ s/:$//;
	    &datePrint( "Need to wait for the following jobs to finish:\n$result");
	    open (SH, ">$outputDirectory/$project/scripts/PeakCallingDependentScript.sh");
	    print SH $header;
	    print SH "#MSUB -W depend=afterok:$result\n";
	    print SH "#MSUB -N CheckingPeakCallerProgress\n";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    print SH "\necho \"Peaking calling jobs $result have finished.\"\n";
	    close SH;
	    &datePrint("Creating dependent job that will only run after peak callers finish.");
	    my $result2 = `msub $outputDirectory/$project/scripts/PeakCallingDependentScript.sh`;
	    $result2 =~ s/\s+//g;
	    my $jobfinished = "no";
	    # Wait until the job is Complete.
	    &datePrint("Waiting for job $result2 to finish. (each . = 30 seconds)");
	    # Check qstat every 30 seconds, adding a "." every time you check.
	    until ($jobfinished eq "Completed"){
		$jobfinished = `checkjob $result2 | grep ^State:`;
		sleep(30);
		print STDERR ".";
		if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
		print STDERR "$jobfinished";
	    }
	    print STDERR "\n";
	    &datePrint("Job $result2 done.  Continuing.");
	}
    }
}

if (($buildDiffPeaks ==1) && ($type eq "chipseq")){
    foreach my $project (keys(%samples)){
	&datePrint("Initiating differential peak analysis for project $project");
	if ($startFromBAM == 0){$bamDirectory = "$outputDirectory\/$project\/bam";}
	if ($chipDescription ne ""){
	    $cmd = "mkdir $outputDirectory\/$project\/analysis";
	    system($cmd);
	    open(CHIP,$chipDescription);
	    my %bamfiles;
	    my %replicates;
	    my %peaksets;
	    my %filename;
	    my %peaksetScript;
	    &datePrint("Using $chipDescription file to match replicates for each sample.");
	    while(<CHIP>){
		if ($_ !~ /^IP/){
		    chomp $_;
		    my ($ip,$input,$peakType,$sample,$peakset) = split(/\,/,$_);
		    $peakType = lc $peakType;
		    my $bedfile = "";
		    if ($peakType eq "narrow"){
			$bedfile = "$outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed";
		    } elsif ($peakType eq "broad"){
			$bedfile = "$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed";
		    }
		    my $samplefile = "$outputDirectory\/$project\/peaks\/$sample.bed";
		    if ($sample ne ""){
			if (exists ($replicates{$sample})){
			    $replicates{$sample} .= ",$bedfile";
			}else {  $replicates{$sample} = $bedfile;}
			my @peaksets = split(/\;/,$peakset);
			foreach my $peakset (@peaksets){			
			    if (exists($peaksets{$peakset})){
				$peaksets{$peakset} .= ",$samplefile";
			    }else {  
				$peaksets{$peakset} = $samplefile;
			    }
			    my $peaksetfile =  "$bamDirectory\/$ip.bam";
			    if (exists($bamfiles{$peakset})){
				$bamfiles{$peakset} .= ",$peaksetfile";
			    } else {
				$bamfiles{$peakset} = "$peaksetfile";
			    }
			    $filename{$bedfile} = $ip;
			    $filename{$samplefile} = $sample;
			    $filename{$peaksetfile} = $ip;
			}		    
		    }
		}
	    }
	    close(CHIP);
	    foreach my $peakset (keys(%peaksets)){
		&datePrint("Set up differential peak analysis script for peakset $peakset.");
		open (SH, ">$outputDirectory/$project/scripts/$peakset\_diffPeaks.sh");
		print SH $header;
		print SH "#MSUB -N $peakset\_diffPeak\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		print SH "module load bedtools/2.17.0\n";
		print SH "module load samtools/1.2\n";
		print SH "module load R/3.2.2\n";
		my @samples = uniq(split(/\,/,$peaksets{$peakset}));
		print SH "\n# Samples for this peakset $peakset are @samples\n";
		foreach my $sample (@samples){
		    $sample = $filename{$sample};
		    my @replicates = split(/\,/,$replicates{$sample});
		    if ($#replicates > 0){
#			&datePrint("Finding reproducible peaks for $sample in $replicates{$sample}");
			print SH "\n# Finding reproducible peaks for $sample\n";
#			print SH "# bedtools multiinter is used, so that replicate numbers greater than 2 can be processed without problem.\n";
#			print SH "# The problem with multiinter is that it takes only the part of the peak that overlaps in every file, and this isn't flexible.  May need to switch to bedtools intersect instead.\n";
#			print SH "bedtools multiinter -i @replicates | \\\n\t";
#			print SH "grep -vw 0 | \\\n\t";
			#  print SH "awk \'{printf \"\%s\\t\%d\\t\%d\\t\%s\_\%d\\t\%d\\n\", \$1,\$2,\$3,\"$sample\",NR,1000}\' \\\n\t";
			my $intersectstring = "bedtools intersect -a $replicates[0] -b $replicates[1] \\\n\t";
			if ($#replicates > 1){
			    # This is supposed to extend the functionality in case of three or more replicates.  It has not been tested.
			    for (my $i=2;$i<=($#replicates+1);$i++){
				$intersectstring .= "| bedtools intersect -a - -b $replicates[$i] \\\n\t";
			    }
			}
			print SH $intersectstring;
			print SH "> $outputDirectory\/$project\/peaks\/$sample.bed\n";
			print SH "# If peaks are close together (within 1 kb), merge them.\n";
			print SH "sort -k1,1 -k2,2n $outputDirectory\/$project\/peaks\/$sample.bed > $outputDirectory\/$project\/peaks\/$sample.sorted.bed\n";
			print SH "bedtools merge -d 1000 -i $outputDirectory\/$project\/peaks\/$sample.sorted.bed > $outputDirectory\/$project\/peaks\/$sample.clustered.bed\n";
			print SH "date\n";
		    } else {
			print SH "\n# Only one replicate for sample $sample, so just creating links.\n";
			print SH "ln -s $replicates{$sample} $outputDirectory\/$project\/peaks\/$sample.bed\n";
		    }
		}
		@samples = split(/\,/,$peaksets{$peakset});
		if ($#samples > 0){
#		    &datePrint("Merging peak files for all samples in the same experimental peak set $peakset: @samples");
		    print SH "\n# Merging peak files for all samples in the same experimental peak set $peakset\n";
		    print SH "cat @samples | \\\n\tbedtools sort -i - | bedtools merge -i - | \\\n\t";
		    print SH "awk \'{printf \"\%s\\t\%d\\t\%d\\t\%s\_\%d\\t\%d\\n\", \$1,\$2,\$3,\"$peakset\",NR,1000}\' \\\n\t";		    
		    print SH "> $outputDirectory\/$project\/peaks\/$peakset.bed\n";
		    print SH "date\n";
		} else {
		    print SH "\n# Only one sample for peakset $peakset, so just creating links.\n";
		    print SH "ln -s $peaksets{$peakset} $outputDirectory\/$project\/peaks\/$peakset.bed\n";
		}
		my @bamfiles = split(/\,/,$bamfiles{$peakset});
		print SH "\n# Calculate coverage for all samples corresponding to the peakset $peakset.\n";
#		&datePrint("Calculate coverage for all samples corresponding to the peakset $peakset.");		
		my @coveragefiles;
		my @coveragefiles2;
		my $tableHeader = "chr\tstart\tstop\tname";
		system("sh $NGSbartom/tools/moduleLoadSamtools.sh");
		foreach my $bamfile (@bamfiles){
		    my $outputfile = "$outputDirectory\/$project\/peaks\/$peakset.".$filename{$bamfile}.".counts.bed";
		    my $outputfile2 = "$outputDirectory\/$project\/peaks\/$peakset.".$filename{$bamfile}.".cpm.bed";
		    print SH "\n# Calculate coverage for $peakset in $bamfile.\n";
		    print SH "bedtools coverage -counts -abam $bamfile -b $outputDirectory\/$project\/peaks\/$peakset.bed | \\\n\t";
		    print SH "awk \'{printf \"\%s\\t\%d\\t\%d\\t\%s\\t\%d\\n\", \$1,\$2,\$3,\$4,\$6}\' | \\\n\t";
		    print SH "sort -k 1,1 -k2,2n > $outputfile\n";
		    print SH "date\n";
		    &datePrint("Checking for number of reads in $bamfile in order to calculate cpm.");
		    my $readCount = `samtools flagstat $bamfile | grep \"+ 0 mapped\" | awk \'{print \$1}\'`;
		    chomp $readCount;
		    print SH "# Num Mapped Reads in $bamfile = $readCount\n";
		    &datePrint("Num Mapped Reads in $bamfile = $readCount");
		    print SH "# Calculate normalized read counts for $peakset in $bamfile.\n";
		    print SH "bedtools coverage -counts -abam $bamfile -b $outputDirectory\/$project\/peaks\/$peakset.bed | \\\n\t";
		    print SH "awk \'{printf \"\%s\\t\%d\\t\%d\\t\%s\\t\%5.3f\\n\", \$1,\$2,\$3,\$4,\(\(\$6*1000000\)/$readCount\)}\' | \\\n\t";
		    print SH "sort -k 1,1 -k2,2n > $outputfile2\n";
		    print SH "date\n";
		    $tableHeader.= "\t$filename{$bamfile}";
		    push (@coveragefiles,$outputfile);
		    push (@coveragefiles2,$outputfile2);
		}
		print SH "\n# Join together counts from all of the bam files for peakset $peakset\n";
		print SH "perl $NGSbartom\/tools\/makePeakCountsTable.pl $outputDirectory\/$project\/peaks\/ $outputDirectory\/$project\/analysis\/ $peakset\n";
		print SH "\n# Generate MDS plot for peakset $peakset\n";
		print SH "Rscript $NGSbartom/tools/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --outputDirectory=$outputDirectory\/$project\/analysis/ --countFile=$outputDirectory\/$project\/analysis\/$peakset.all.counts.txt --numCores=$numProcessors --runMDS=1\n";
		&datePrint("Looking for $outputDirectory\/$project\/$peakset.comparisons.csv");
		if (-e "$outputDirectory\/$project\/$peakset.comparisons.csv"){
		    print SH "\n# Find Differential Peaks for peakset $peakset\n";
		    print SH "Rscript $NGSbartom/tools/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --outputDirectory=$outputDirectory\/$project\/analysis/ --countFile=$outputDirectory\/$project\/analysis\/$peakset.all.counts.txt --numCores=$numProcessors --comparisonFile=$outputDirectory\/$project\/$peakset.comparisons.csv --runMDS=0\n";

		    ## For each comparison, go through and output Bed file of significant peaks
		} else {
		    &datePrint("To find differentially expressed peaks in peakset $peakset, create a comparisons file at $outputDirectory\/$project\/$peakset.comparisons.csv")
		}
		print SH "date\n";
		close(SH);
		if ($runDiffPeaks == 1){
		    &datePrint("Starting job $outputDirectory/$project/scripts/$peakset\_diffPeaks.sh");
		    `msub $outputDirectory/$project/scripts/$peakset\_diffPeaks.sh`;
		}
	    }
	} else {
	    die "ERR:  ChipDescription file must be provided which specifies which samples are replicates and which are in the same peakset.\n";
	}   
    }
    # if ($runDiffPeaks ==1){
    # 	# Submit the diffPeak jobs, saving the job id.
    # 	## This will finish each project (TANGO) before starting the next.  Is this what we want?
    # 	foreach my $project (keys(%samples)){
    # 	    &datePrint ("Starting diff peak scripts.");
    # 	    my $result = `find $outputDirectory/*/scripts/ -iname \"*diffPeaks.sh\" -exec msub {} ./ \\\;`;
    # 	    $result =~ s/\s+/\:/g;
    # 	    $result =~ s/^://;
    # 	    $result =~ s/:$//;
    # 	    &datePrint( "Need to wait for the following jobs to finish:\n$result");
    # 	    open (SH, ">$outputDirectory/$project/scripts/diffPeakDependentScript.sh");
    # 	    print SH $header;
    # 	    print SH "#MSUB -W depend=afterok:$result\n";
    # 	    print SH "#MSUB -N CheckingDiffPeakProgress\n";
    # 	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
    # 	    print SH "\necho \"Peaking calling jobs $result have finished.\"\n";
    # 	    close SH;
    # 	    &datePrint("Creating dependent job that will only run after diff peak scripts finish.");
    # 	    my $result2 = `msub $outputDirectory/$project/scripts/diffPeakDependentScript.sh`;
    # 	    $result2 =~ s/\s+//g;
    # 	    my $jobfinished = "no";
    # 	    # Wait until the job is Complete.
    # 	    &datePrint("Waiting for job $result2 to finish. (each . = 30 seconds)");
    # 	    # Check qstat every 30 seconds, adding a "." every time you check.
    # 	    until ($jobfinished eq "Completed"){
    # 		$jobfinished = `checkjob $result2 | grep ^State:`;
    # 		sleep(30);
    # 		print STDERR ".";
    # 		if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
    # 		print STDERR "$jobfinished";
    # 	    }
    # 	    print STDERR "\n";
    # 	    &datePrint("Job $result2 done.  Continuing.");
    #}
    #}
}


&datePrint("Finished.");

sub datePrint{
    my $printString = $_[0];
    my $date = `date --rfc-3339='ns'`;
    if ($date =~ /(\d+\-\d+\-\d+)\s(\d+\:\d+\:\d+)\./){
	$date = "[$1 $2]";
    }
#    print STDERR "Date = $date\n";
    print STDERR "$date\t$printString\n";
}
