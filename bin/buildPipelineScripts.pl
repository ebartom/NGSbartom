#!/software/activeperl/5.16/bin/perl -w
use Getopt::Long qw(GetOptions);
use List::Util qw(max);
use List::MoreUtils qw(uniq);
use File::Basename;
use strict;
use utf8;
use warnings;
no warnings 'uninitialized';
#use Config::Abstract::Ini;

unless (@ARGV) {
    print "\nUsage: buildPipelineScripts.pl\n\n-config\t\t<configFile>\n-o\t\t<outputDirectory>\n-bs\t\t<baseSpaceDirectory>\n-f\t\t<fastqDirectory>\n";
    print "-bam\t\t<bamDirectory>\n-c\t\tcomparisons.csv file>\n-4C\t\t<4C description file>\n-chip\t\t<ChIP Description file>\n";
    print "-p\t\t<numProcessors>\n-m\t\t<multiMap (1 = allow multi map, 0 = not)>\n";
    print "-a\t\t<aligner>\n-g\t\t<assembly/genome>\n-w\t\t<walltime>\n-t\t\t<RNA|chipseq|4C>\n-account\t<accountName>\n";
    print "-node\t\t<nodeName(eg qnode4144)>\n-scientist\t<initials>\n-s3\t\t<amazon s3 path (if non-standard)>\n-buildBcl2fq\t\t<1|0>\n";
    print "-runBcl2fq\t\t<1|0>\n-runTrim\t\t<1|0>\n-buildAlign\t\t<1|0>\n-runAlign\t\t<1|0>\n-makeTracks\t\t<1|0>\n";
    print "-uploadASHtracks\t<1|0>\n-uploadBAM\t\t<1|0>\n-buildEdgeR\t\t<1|0>\n-runEdgeR\t\t<1|0>\n-buildPeakCaller\t<1|0>\n-runPeakCaller\t\t<1|0>\n";
    print "-buildDiffPeaks\t\t<1|0>\n-runDiffPeaks\t\t<1|0>\n-build4C\t\t<1|0>\n-run4C\t\t\t<1|0>\n-runRNAstats\t\t<1|0>\n";
}

# Set up the environment for analysis.
my $NGSbartom="/projects/p20742/";

my $configFile = "";
my $distToTSS = 2000;
my $upstream = 5000;
my $downstream = 5000;
my $trimString = "TRAILING:30 MINLEN:20";
my $tophatReadMismatch = 2; # 2 is Tophat default.
my $tophatReadEditDist = 2; # 2 is Tophat default.
my $tophatMultimap = 20; # 20 is Tophat default.
my $runPairedEnd = 0;
my $outputDirectory = "";
my $baseSpaceDirectory = "";
my $bamDirectory = "";
my $fastqDirectory = "";
my $sampleSheet = "";
my $comparisons = "";
my $type = "";
my $scheduler = "SLURM";
my $numProcessors = 4;
my $memory = 50000;
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
my $uploadLutracks = 0;
my $uploadBAM = 0;
my $buildEdgeR = 0;
my $runTrim = 1;
my $buildBcl2fq = 0;
my $buildAlign = 0;
my $buildSampleCheck = 0;
my $runSampleCheck = 0;
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
my $exomeDescription = "";
my $htseq = 1;
my $bedtools = 0;
my $ngsplot = 0;
my $granges = 0;
my $rsem = 0;
my $genomeBAM = 1;
my $rgString = "";
my $runRNAstats = 0;
my $buildNGSplot = 0;
my $ngsFCcomparison = "";

# S3path is not working properly right now for non-standard s3path values.  This needs fixing.

## NB: Below there's an option to specify a read group string, but that would
## set the read group to be the same for all samples, which doesn't make sense.
## So for now, zero-ing out readgroup when we start iterating through samples.


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
	   'exomePairingFile|pairs=s' => \$exomeDescription,
	   'type|t=s' => \$type,
	   'aligner|a=s' => \$aligner,
	   'trimString|ts=s' => \$trimString,
	   'assembly|g=s' => \$assembly,
	   'tango|id=s' => \$id,
	   'scientist|s=s' => \$scientist,
	   'stranded|str=i' => \$stranded,
	   'processors|p=i' => \$numProcessors,
	   'memory|mem=i' => \$memory,
	   'scheduler|sch=s' => \$scheduler,
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
	   'uploadLutracks|ult=i' => \$uploadLutracks,
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
	   'buildSampleCheck|bsc=i' => \$buildSampleCheck,
	   'runSampleCheck|rsc=i' => \$runSampleCheck,
	   'build4C|b4=i' => \$build4C,
	   'run4C|r4=i' => \$run4C,
	   'runBcl2fq|rb=i' => \$runBcl2fq,
	   'htseq|h=i' => \$htseq,
	   'bedtools|bt=i' => \$bedtools,
	   'ngsplot|np=i' => \$ngsplot,
	   'granges|ash=i' => \$granges,
	   'rsem=i' => \$rsem,
	   'genomeBAM|gb=i' => \$genomeBAM,
	   'runRNAstats|rrs=i' => \$runRNAstats,
	   'buildNGSplot|ngsp=i' => \$buildNGSplot,
	   'ngsFCcomparison|ngspc' => \$ngsFCcomparison
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
if (($uploadLutracks == 1) && ($uploadASHtracks == 1)){
    print STDERR "ERROR: Set to upload to both Shilatifard S3 account and Lu S3 account.  Since Shilatifard is default and Lu is manually set, assuming only Lu is desired and resetting flags accordingly.\n";
    $uploadASHtracks = 0;
}
if (($s3path eq "") && ($uploadASHtracks == 1)){ $s3path = "ash-tracks/TANGO/$scientist";}
if (($s3path eq "") && ($uploadPulmtracks == 1)){ $s3path = "m-328-data/$scientist";}
if (($s3path eq "") && ($uploadLutracks == 1)){ $s3path = "lwwlab/";}
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
    if ($type eq "proseq"){ $aligner = "bowtie";}
    if ($type eq "4su"){ $aligner = "bowtie";}
    if ($type eq "4C"){ $aligner = "bowtie";}
} else { $aligner = lc $aligner;}

if (($buildSampleCheck == 0) && ($buildAlign == 1) && (($assembly eq "hg19") || ($assembly eq "hg38") || ($assembly eq "hg38.mp"))){
    print STDERR "If your assembly is human and you are doing an alignment, turning on sample check by default.\n";
    $buildSampleCheck = 1;
}
if (($buildSampleCheck == 1) && 
    (($assembly eq "mm9") || ($assembly eq "mm10") || ($assembly eq "dm3") || ($assembly eq "sacCer3") || ($assembly eq "rn6") || ($assembly eq "grcm38"))){
    print STDERR "ERROR: buildSampleCheck currently only works for human samples. Turning it off.\n";
    $buildSampleCheck = 0;
}   

# You can build scripts without running them, but not the reverse.
if (($runBcl2fq == 1) && ($buildBcl2fq == 0)){
    $buildBcl2fq = 1;
}
if (($runAlign == 1) && ($buildAlign == 0)){
    $buildAlign = 1;
}
if (($runSampleCheck == 1) && ($buildSampleCheck == 0)){
    $buildSampleCheck = 1;
}
if (($runGenotyping == 1) && ($buildGenotyping == 0)){
    $buildGenotyping = 1;
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

$scheduler = uc($scheduler);
print STDERR "Setting up scripts for $scheduler scheduler.\n";
print STDERR "Account: $account\nQueue/Partition: $queue\n";

my $header = "#!/bin/bash\n";
if ($scheduler eq "MOAB"){
    # Define a header for the shell scripts.
    $header .= "#MSUB -A $account\n";
    if ($account eq "b1042"){
	$header .= "#MSUB -q $queue\n";
    } elsif ($account eq "e30258"){
	$header .= "#MSUB -q $queue\n";
    } elsif ($account eq "b1025"){
	if ($queue eq ""){ $queue = "buyin";}
	$header .= "#MSUB -q $queue\n";
    }
    $header .= "#MSUB -l walltime=$walltime\n";
    $header .= "#MSUB -m a\n"; # only email user if job aborts
    $header .= "#MSUB -j oe\n";
    # This umask parameter doesn't seem to make sense so commenting it out for now.
    #    $header .= "#MOAB -W umask=0113\n";
    if ($node ne ""){
	$header .= "#MSUB -l nodes=$node\n";
    }
} elsif ($scheduler eq "SLURM"){
    $header .= "#SBATCH -A $account\n";
    if ($account eq "b1042"){
	$header .= "#SBATCH -p $queue\n";
    } elsif ($account eq "e30258"){
	$header .= "#SBATCH -p $queue\n";
    } elsif ($account eq "e30682"){
	$header .= "#SBATCH -p $queue\n";
    } elsif ($account eq "b1025"){
	if ($queue eq ""){ $queue = "buyin";}
	$header .= "#SBATCH -p $queue\n";
    }
    $header .= "#SBATCH -t $walltime\n";
    $header .= "#SBATCH -m a\n"; # only email user if job aborts
    $header .= "#SBATCH --mem=$memory\n";
    # Run the scripts from $outputDirectory/metadata, which will put the stdout / stderr files there for debugging.
    # First check for metadata directory existence, and create if absent
    if (!(-e "$outputDirectory\/metadata")){
      `mkdir $outputDirectory\/metadata`;
    }
    $header .= "#SBATCH --chdir=$outputDirectory/metadata/\n";
    $header .= "#SBATCH -o \"\%x.o\%j\"\n";
    #    $header .= "#MOAB -W umask=0113\n"; #How do I do this in SLURM?
    # Maybe we just shouldn't be doing this.  Leaving it out for now.
    if ($node ne ""){
	$header .= "#SBATCH --nodelist=$node\n";
    }
} else {
    print STDERR "ERR: Unrecognized scheduler $scheduler.\n";
}


print STDERR "=====================================\n";
print STDERR "Analysis should be initiated either with a base space directory, or a fastq directory, or a bam directory.  Below you can see what was specified in this run.\n";
print STDERR "BaseSpaceDirectory: $baseSpaceDirectory\n";
print STDERR "Fastq Directory: $fastqDirectory\n";
print STDERR "Bam Directory: $bamDirectory\n";
print STDERR "Output Directory: $outputDirectory\n";
print STDERR "Type of Analysis: $type\n";
print STDERR "=====================================\nAnalysis plan:\n";
if ($runPairedEnd){ print STDERR "All fastq files are assumed to be paired, with pairings indicated by _R1, _R2 or _R1, _R3.\n";}
if ($buildBcl2fq){ print STDERR "Will build scripts for de-multiplexing with Bcl2fq.\n";}
if ($runBcl2fq){ print STDERR "Will run scripts for de-multiplexing with Bcl2fq.\n";}
if ($runAlign && $runTrim) { print STDERR "Will trim fastq files according to trailing quality scores with Trimmomatic ($trimString).\n";}
if ($buildAlign){ print STDERR "Will build scripts for aligning fastq reads with $aligner.\n";}
if ($runAlign){ print STDERR "Will run scripts for aligning fastq reads with $aligner.\n";}
if ($buildSampleCheck){ print STDERR "Will build scripts for checking sample identity for $assembly\n";}
if ($runSampleCheck){ print STDERR "Will run scripts for checking sample identity for $assembly.\n";}
if ($buildPeakCaller){ print STDERR "Will build scripts for calling peaks according to the experimental plan in $chipDescription.\n";}
if ($runPeakCaller){ print STDERR "Will run scripts for calling peaks.\n";}
if ($makeTracks){ print STDERR "Will make tracks appropriate for the UCSC genome browser.\n";
		  if ($uploadASHtracks){ print STDERR "Will upload tracks to Shilatifard account on Amazon S3 (if user has correct credentials).\n";}
		  if ($uploadPulmtracks){ print STDERR "Will upload tracks to Pulmonology account on Amazon S3 (if user has correct credentials).\n";}
		  if ($uploadLutracks){ print STDERR "Will upload tracks to Lu Wang's Lab account on Amazon S3 (if user has correct credentials).\n";}
		  if ($uploadBAM){ print STDERR "Will upload BAM files to Shilatifard account on Amazon S3 (if user has correct credentials).\n";}
}
if ($runRNAstats) { print STDERR "Will launch a script to assess RNAseq quality based on BAM files.\n";}
if ($buildEdgeR){ print STDERR "Will build scripts for finding differentially expressed genes from RNAseq data.\n";}
if ($runEdgeR){ print STDERR "Will run scripts for finding differentially expressed genes from RNAseq data.\n";}
if ($buildDiffPeaks) { print STDERR "Will build scripts to find differentially expressed peaks based on $chipDescription and $comparisons.\n";}
if ($runDiffPeaks) { print STDERR "Will run scripts to find differentially expressed peaks.\n";}
if ($build4C){ print STDERR "Will build scripts for analyzing 4C data, starting with mock-demultiplexed fastq files.\n";}
if ($run4C){ print STDERR "Will run scripts for analyzing 4C data.\n";}
if ($buildGenotyping){ print STDERR "Will build scripts for genotyping samples.\n";}
if ($runGenotyping){ print STDERR "Will run scripts for genotyping samples.\n";}
print STDERR "=====================================\n";

# Define references.
my (%bowtieIndex,%starIndex,%txIndex,%txdbfile,%bwaIndex,%gff,%exonbed,%rsemTx,%genebed,%samplegenebed);
my (%gatkRef,%knownSNPsites,%knownIndelsites,%cosmicSNPsites);

$bowtieIndex{"hg38"} = "$NGSbartom/anno/bowtie_indexes/hg38";
$bwaIndex{"hg38"} = "$NGSbartom/anno/bwa_indexes/hg38.fa";
$starIndex{"hg38"} = "$NGSbartom/anno/STAR_indexes/hg38/";
$txIndex{"hg38"} ="$NGSbartom/anno/tophat_tx/hg38.Ens_78.remap";
$txdbfile{"hg38"} = "$NGSbartom/anno/Txdb/hsapiens_gene_ensembl_Ens78.txdb";
$exonbed{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.exons.bed";
$gff{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.cuff.gtf";
$genebed{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.cuff.bed"; # Created with gtf2bed tool
$samplegenebed{"hg38"} = "$NGSbartom/anno/Ens/hg38.Ens_78/hg38.Ens_78.r1k.cuff.bed"; # 1000 genes picked at random from gene bed.
$rsemTx{"hg38"} = "$NGSbartom/anno/rsemTx/hg38.Ens_78";
$gatkRef{"hg38"} = "$NGSbartom/anno/picardDict/hg38.fa";
$knownSNPsites{"hg38"} = "$NGSbartom/anno/picardDict/1000G_phase1.snps.high_confidence.hg38.vcf";
$cosmicSNPsites{"hg38"} = "$NGSbartom/anno/Cosmic/CosmicCodingMuts.vcf.gz";
$knownIndelsites{"hg38"} = "$NGSbartom/anno/picardDict/Mills_and_1000G_gold_standard.indels.hg38.vcf";

$bowtieIndex{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome";
$bwaIndex{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa";
$starIndex{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Sequence/STARindex/";
$txIndex{"hg38.mp"} ="$NGSbartom/anno/tophat_tx/hg38.Ens_78.remap";
$txdbfile{"hg38.mp"} = "$NGSbartom/anno/Txdb/UCSC.hg38.mp.txdb";
$exonbed{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf.bed";
$genebed{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf.bed";
$samplegenebed{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.r1k.gtf.bed"; # 1000 genes picked at random from gene bed.
$gff{"hg38.mp"} = "$NGSbartom/anno/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf";
$rsemTx{"hg38.mp"} = "NEEDS TO BE UPDATED $NGSbartom/anno/rsemTx/hg38.Ens_78";
$gatkRef{"hg38.mp"} = "$NGSbartom/anno/picardDict/hg38.fa";
$knownSNPsites{"hg38.mp"} = "$NGSbartom/anno/picardDict/1000G_phase1.snps.high_confidence.hg38.vcf";
$cosmicSNPsites{"hg38"} = "$NGSbartom/anno/Cosmic/CosmicCodingMuts.vcf.gz";
$knownIndelsites{"hg38.mp"} = "$NGSbartom/anno/picardDict/Mills_and_1000G_gold_standard.indels.hg38.vcf";

$bowtieIndex{"hg19"} = "$NGSbartom/anno/bowtie_indexes/hg19";
$bwaIndex{"hg19"} = "$NGSbartom/anno/bwa_indexes/hg19.fa";
$starIndex{"hg19"} = "$NGSbartom/anno/STAR_indexes/hg19_ens75/";
$txIndex{"hg19"} ="$NGSbartom/anno/tophat_tx/hg19.Ens_75.remap";
$txdbfile{"hg19"} = "$NGSbartom/anno/Txdb/hsapiens_gene_ensembl_Ens75.txdb";
$exonbed{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_75/hg19.Ens_75.exons.bed";
$genebed{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_75/hg19.Ens_75.cuff.bed";
$samplegenebed{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_75/hg19.Ens_75.r1k.cuff.bed"; # 1000 genes picked at random from gene bed.
$gff{"hg19"} = "$NGSbartom/anno/Ens/hg19.Ens_75/hg19.Ens_75.cuff.gtf";
$rsemTx{"hg19"} = "$NGSbartom/anno/rsemTx/hg19.Ens_75";
$gatkRef{"hg19"} = "$NGSbartom/anno/picardDict/hg19.fa";
$knownSNPsites{"hg19"} = "$NGSbartom/anno/picardDict/1000G_phase1.snps.high_confidence.hg19.sites.noContigs.vcf";
$knownIndelsites{"hg19"} = "$NGSbartom/anno/picardDict/1000G_phase1.indels.hg19.sites.noContigs.vcf";

$bowtieIndex{"dm3"} = "$NGSbartom/anno/bowtie_indexes/dm3";
$bwaIndex{"dm3"} = "$NGSbartom/anno/bwa_indexes/dm3.fa";
$starIndex{"dm3"} = "$NGSbartom/anno/STAR_indexes/dm3/";
$txIndex{"dm3"} = "$NGSbartom/anno/tophat_tx/dm3.Ens_74.cuff";
$txdbfile{"dm3"} = "$NGSbartom/anno/Txdb/dmelanogaster_gene_ensembl_Ens74.txdb";
$exonbed{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.exons.bed";
$genebed{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.cuff.bed";
$samplegenebed{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.r1k.cuff.bed";
$gff{"dm3"} = "$NGSbartom/anno/Ens/dm3.Ens_74/dm3.Ens_74.cuff.gtf";
$rsemTx{"dm3"} = "$NGSbartom/anno/rsemTx/dm3.Ens_74";

$bowtieIndex{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Sequence/Bowtie2Index/genome";
$bwaIndex{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Sequence/BWAIndex/genome.fa";
$txIndex{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.txIndex";
$txdbfile{"grcm38"} = "$NGSbartom/anno/Txdb/mmusculus_gene_ensembl_Ens87.txdb";
$exonbed{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/exons.bed";
$gff{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf";
$rsemTx{"grcm38"} = "$NGSbartom/anno/rsemTx/mouse_GRCm38";
$gatkRef{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa";
$knownSNPsites{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Annotation/Variation/Mus_musculus.reordered2.vcf";
$knownIndelsites{"grcm38"} = "$NGSbartom/anno/Mus_musculus/Ensembl/GRCm38/Annotation/Variation/Mus_musculus.reordered2.vcf";

$bowtieIndex{"rn6"} = "$NGSbartom/anno/bowtie_indexes/rn6";
$bwaIndex{"rn6"} = "$NGSbartom/anno/bwa_indexes/rn6.fa";
$starIndex{"rn6"} = "$NGSbartom/anno/STAR_indexes/rn6/";
$txIndex{"rn6"} = "$NGSbartom/anno/tophat_tx/rn6.UCSC";
$txdbfile{"rn6"} = "$NGSbartom/anno/Txdb/rnorvegicus_gene_ensembl_Ens84.txdb";
$exonbed{"rn6"} = "$NGSbartom/anno/Ens/rn6.Ens_84/rn6.Ens_84.exons.bed";
$genebed{"rn6"} = "$NGSbartom/anno/Ens/rn6.Ens_84/rn6.Ens_84.cuff.bed";
$samplegenebed{"rn6"} = "$NGSbartom/anno/Ens/rn6.Ens_84/rn6.Ens_84.r1k.cuff.bed";
#$gff{"rn6"} = "$NGSbartom/anno/Ens/rn6.Ens_84/rn6.Ens_84.cuff.gtf";
#$gff{"rn6"} = "$NGSbartom/anno/Ens/rn6.Ens_84/Rattus_norvegicus.Rnor_6.0.84.gtf";
$gff{"rn6"} = "/projects/p20742/anno/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf";
$rsemTx{"rn6"} = "$NGSbartom/anno/rsemTx/rn6.Ens_84";

$bowtieIndex{"mm10"} = "$NGSbartom/anno/bowtie_indexes/mm10";
$bwaIndex{"mm10"} = "$NGSbartom/anno/bwa_indexes/mm10.fa";
$starIndex{"mm10"} = "$NGSbartom/anno/STAR_indexes/mm10/";
$txIndex{"mm10"} = "$NGSbartom/anno/tophat_tx/mm10.Ens_78.cuff";
$txdbfile{"mm10"} = "$NGSbartom/anno/Txdb/mmusculus_gene_ensembl_Ens78.txdb";
$exonbed{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.exons.bed";
$genebed{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.bed";
$samplegenebed{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.r1k.cuff.bed";
$gff{"mm10"} = "$NGSbartom/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf";
$rsemTx{"mm10"} = "$NGSbartom/anno/rsemTx/mm10.Ens_78";

$bowtieIndex{"mm9"} = "$NGSbartom/anno/bowtie_indexes/mm9";
$bwaIndex{"mm9"} = "$NGSbartom/anno/bwa_indexes/mm9.fa";
$starIndex{"mm9"} = "$NGSbartom/anno/STAR_indexes/mm9/";
$txIndex{"mm9"} = "$NGSbartom/anno/tophat_tx/mm9.Ens_67.remap";
$txdbfile{"mm9"} = "$NGSbartom/anno/Txdb/mmusculus_gene_ensembl_Ens67.txdb";
$exonbed{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.exons.bed";
$genebed{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.cuff.bed";
$samplegenebed{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.r1k.cuff.bed";
$gff{"mm9"} = "$NGSbartom/anno/Ens/mm9.Ens_67/mm9.Ens_67.cuff.gtf";
$rsemTx{"mm9"} = "$NGSbartom/anno/rsemTx/mm9.Ens_67";

$bowtieIndex{"sacCer3"} = "$NGSbartom/anno/bowtie_indexes/sacCer3";
$bwaIndex{"sacCer3"} = "$NGSbartom/anno/bwa_indexes/sacCer3.fa";
$starIndex{"sacCer3"} = "$NGSbartom/anno/STAR_indexes/sacCer3/";
$txIndex{"sacCer3"} = "$NGSbartom/anno/tophat_tx/sacCer3.Ens_72.remap";
#$txdbfile{"sacCer3"} = "$NGSbartom/anno/Txdb/scerevisiae_gene_ensembl_Ens72.txdb";
#$txIndex{"sacCer3"} = "$NGSbartom/anno/tophat_tx/sacCer3.Ens_78.remap";
$txdbfile{"sacCer3"} = "$NGSbartom/anno/Txdb/scerevisiae_gene_ensembl_Ens78.txdb";
$exonbed{"sacCer3"} = "$NGSbartom/anno/Ens/sacCer3.Ens_78/sacCer3.Ens_78.exons.bed";
$genebed{"sacCer3"} = "$NGSbartom/anno/Ens/sacCer3.Ens_78/sacCer3.Ens_78.cuff.bed";
$samplegenebed{"sacCer3"} = "$NGSbartom/anno/Ens/sacCer3.Ens_78/sacCer3.Ens_78.r1k.cuff.bed";
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

my @now = localtime();
my $timestamp = sprintf("%04d.%02d.%02d.%02d%02d.%02d", $now[5]+1900, $now[4]+1, $now[3],$now[2],$now[1],$now[0]);
my $timestamp2 = sprintf("%04d-%02d-%02d", $now[5]+1900, $now[4]+1, $now[3]);
#print STDERR "Timestamp $timestamp\n";

# If buildBcl2fq == 1, then create a shell script for Bcl2fq (if runBcl2fq == 1, then submit the job and wait for it to finish).
if ($buildBcl2fq == 1){
    &datePrint("Creating shell script for Bcl2fq");
    my $shScript = "$baseSpaceDirectory\/runBcl2fq.sh";
    my $bclFqProcessors = max($numProcessors,8);
    if (!(-e "$outputDirectory\/metadata")){
	`mkdir $outputDirectory\/metadata`;
    }
    open(VER,">$outputDirectory\/metadata\/Ceto.run.$type.$timestamp.txt");
    open(SH,">$shScript");
    print SH $header;
    if ($scheduler eq "MOAB"){
	print SH "#MSUB -l nodes=1:ppn=$bclFqProcessors\n";
	print SH "#MSUB -N bcl2fastq\n";
    } elsif ($scheduler eq "SLURM"){
	print SH "#SBATCH --nodes=1\n";
	print SH "#SBATCH -n $bclFqProcessors\n";
	print SH "#SBATCH --job-name=bcl2fastq\n";
    }
    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
    print SH "\nmodule load bcl2fastq/2.17.1.14\n";
    print VER "module load bcl2fastq/2.17.1.14\n";
    print SH "bcl2fastq -R $baseSpaceDirectory -r $numProcessors -d $numProcessors -p $numProcessors -w $numProcessors\n";
    # If runBcl2fq == 1, then run the shell script (only works if buildBcl2fq == 1)
    if ($runBcl2fq == 1){
	&datePrint("Running Bcl2fq job and waiting for it to finish.");
	# Submit the job, saving the job id.
	my $result = "";
	if ($scheduler eq "MOAB"){
	    $result = `msub $shScript`;
	    $result =~ s/\s//g;
	} elsif ($scheduler eq "SLURM"){
	    $result = `sbatch $shScript`;
	    $result =~ s/Submitted batch job //g;
#	    print STDERR "Predicted SLURM job ID:  $result\n";
	}
	my $jobfinished = "no";
	# Wait until the job is no longer running.
	&datePrint("Waiting for job $result to finish.");
	waitForJob($result);
	print STDERR "\n";
	&datePrint("Job $result done.  Continuing.");
    }
    close SH;
    $cmd = "cp $baseSpaceDirectory\/runBcl2fq.sh $outputDirectory\/\n";
    $cmd .= "cp $baseSpaceDirectory\/SampleSheet.csv $outputDirectory\/metadata\/\n";
    system($cmd);
}
my %scientists;

my(%fastqs,%reference,%samples,%sampleIDs,$sample_project,$sample_plate);
$sample_plate = $scientist;
if (($sampleSheet ne "")){
    my($sampleIDindex,$sampleNameIndex,$assemblyIndex,$plateIndex,$descriptionIndex,$projectIndex);
# Read in the sample sheet
    &datePrint("Reading in $sampleSheet and printing Sample_Report");
    open (IN,$sampleSheet);
    my $flag="header";
    my($sample_ID,$sample_name,$assembly,$I7_Index_ID,$index,$description,$index2,$I5_Index_ID,$stuff);
    my $sampleNum = 0;
    if (!(-e "$outputDirectory\/metadata\/SampleSheet.csv")){
	my $cmd = "cp $baseSpaceDirectory\/SampleSheet.csv $outputDirectory\/metadata\/";
	system($cmd);
    }
    # Create output file for Sample_Report.
    open(OUT,">$outputDirectory/metadata/Sample_Report.csv");
    # Print labels to output file.
    print OUT "Fastq,Sample_Name,Assembly\n";
    # Create run file to print metadata to.
    open(VER,">$outputDirectory\/metadata\/Ceto.run.$type.$timestamp.txt");
    while(<IN>){
	chomp $_;
	if ($_ =~ /^Sample_ID/){
	    # Check each column header to figure out which to use for what
	    my @headers = split(/\,/,$_);
	    for (my $i=0;$i<=$#headers;$i++){
		if ($headers[$i] eq "Sample_ID"){
		    $sampleIDindex = $i;
		    print STDERR "Sample ID index: $sampleIDindex\n";
		} elsif ($headers[$i] eq "Sample_Name"){
		    $sampleNameIndex = $i;
		    print STDERR "Sample Name index: $sampleNameIndex\n";
		} elsif ($headers[$i] eq "assembly" || $headers[$i] eq "organism"){
		    $assemblyIndex = $i;
		} elsif ($headers[$i] eq "Sample_Plate"){
		    $plateIndex = $i; # We actually use this for scientist.
		} elsif ($headers[$i] eq "Description"){
		    $descriptionIndex = $i;
		} elsif ($headers[$i] eq "Sample_Project"){
		    $projectIndex = $i;
		    print STDERR "Sample Project index: $projectIndex\n";
		}
	    }
	}
	if (($flag eq "data") && ($_ !~ /^Sample_ID/) && ($_ !~ /SampleID/)){
	    # Read in the metadata from the sample sheet.
	    my @data = split(/\,/,$_);
	    #	    ($sample_ID,$sample_name,$sample_plate,$assembly,$I7_Index_ID,$index,$sample_project,$description,$stuff) = split(/\,/,$_);
	    $sample_ID = $data[$sampleIDindex];
	    $sample_name = $data[$sampleNameIndex];
	    $assembly = $data[$assemblyIndex];
	    $sample_plate = $data[$plateIndex];
	    $description = $data[$descriptionIndex];
	    print VER "PROJ From Sample sheet this experiment is $description\n";
	    $sample_project = $data[$projectIndex];
	    $sampleNum++;
	    $sampleIDs{$sample_name} = $sample_ID;
	    # Store the assembly for the sample and project.
	    $reference{$sample_name}=$assembly;
	    $reference{$sample_project} = $assembly;
	    $scientists{$sample_project} = $sample_plate;
	    print VER "PROJ $sample_project\n";
	    print STDERR "Project: $sample_project\n";
	    &datePrint("Found the following sample:");
	    print STDERR "$sample_name\tREF:$reference{$sample_name}\t$bowtieIndex{$assembly}\n";
	    if ($type ne "4C"){
		my $fastq = "";
		my $fastq2 = "";
		# For Nova-seq, there are only two lanes, L001 and L002
		foreach my $lane ("L001","L002","L003","L004"){
		    # Build fastq file names.
		    $fastq = "$sample_name\_S$sampleNum\_$lane\_R1\_001.fastq.gz";
		    if ($runPairedEnd == 1){
			$fastq2 = "$sample_name\_S$sampleNum\_$lane\_R2\_001.fastq.gz";
		    }
		    #		    print STDERR "Looking for $fastq\n";
		    print STDERR "Looking for $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq or $outputDirectory\/$sample_project\/fastq\/$fastq\n";
		    # Check if fastq exists in new subdirectory or old.
		    if ((-e "$baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq") || (-e "$outputDirectory\/$sample_project\/fastq\/$fastq")){
			if ((-e "$baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq") ){
			    if (-e "$outputDirectory\/$sample_project\/fastq\/$fastq") {
				#				print STDERR "$fastq exists in both BaseSpace directory and output directory\n";
			    } else {
				print STDERR "$fastq exists in BaseSpace directory but not output directory; moving it over.\n";
				print STDERR "Deleting sample directory in basespace directory.\n";
				my $cmd = "mv $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq $outputDirectory\/$sample_project\/fastq\/\n";
				if ($runPairedEnd == 1){
				    my $cmd = "mv $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\/$fastq2 $outputDirectory\/$sample_project\/fastq\/\n";
				}
				$cmd .= "rmdir $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/$sample_project\/$sample_ID\n";
				system($cmd);
			    }
			}
			if (-e "$outputDirectory\/$sample_project\/fastq\/$fastq") {
			    # Fastq file is in output directory.  All is good.
			} else {
			    # If you can't find Fastq file anywhere, then die.
			    # This would indicate an error in the path or filename, or that the bcl2fq job failed.  It could also indicate that the index was wrong in the sample sheet, and that no reads were assigned to a specific sample.
			    print STDERR "ERR:  Cannot find fastq file! This could be because NovaSeq has fewer lanes than NextSeq.\n";
			}
		    }
		    # Write the sample report file.
		    print OUT "$fastq,$sample_name,$assembly\n";
		    # Create a hash of all read1 fastq files for a given sample_name.
		    if (!exists($fastqs{$sample_name})){
			$fastqs{$sample_name} = "$outputDirectory\/$sample_project\/fastq\/$fastq";
		    }else {$fastqs{$sample_name} .= ",$outputDirectory\/$sample_project\/fastq\/$fastq";}
		}
	    }
	    # Create a hash of all samples in a given sample project (TANGO/MOLNG)
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
    my ($sample_name,$project_name,@fastqlist,@bamlist,$fastqlist,$bamlist);
    if (($assembly eq "") || (($fastqDirectory eq "")&& ($startFromBAM == 0))){
	die "ERR:  If SampleSheet is not specified, assembly and fastqDirectory or bamDirectory must be\n";
    } elsif ($fastqDirectory ne "") {
	&datePrint("Looking for Fastq files in $fastqDirectory.");
	$fastqlist = "";
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
	@fastqlist = split(/\s+/,$fastqlist);
	&datePrint("Found @fastqlist");
#	my $project_name = "thisProject";
	if ($fastqDirectory =~ /([\w\-\_\.]+)\/fastq\/?$/){
	    $project_name = $1;
	} elsif ($fastqDirectory =~ /([\w\-\_\.]+)\/?$/){
	    $project_name = $1;
	    $project_name =~ s/.fastqs//g;
	    $project_name =~ s/.fastq//g;
	    $project_name =~ s/.seqfiles//g;
	}
	&datePrint("Project name is $project_name");
    } elsif ($bamDirectory ne "") {
	&datePrint("Looking for BAM files in $bamDirectory.");
	$bamlist = `ls $bamDirectory\/*.bam`;
	@bamlist = split(/\s+/,$bamlist);
	&datePrint("Found @bamlist");
	$project_name = basename($bamDirectory);
	if ($bamDirectory =~ /([\w\-\_\.]+)\/bam\/?$/){
	    $project_name = $1;
	}
	$project_name =~ s/.bams//g;
	$project_name =~ s/.bam//g;
	$project_name =~ s/.bamfiles//g;
	&datePrint("Project name is $project_name");
    }
    if (($fastqDirectory ne "") || ($bamDirectory ne "")){
	$reference{$project_name}=$assembly;
	my @now = localtime();
	my $timestamp = sprintf("%04d.%02d.%02d.%02d%02d.%02d", $now[5]+1900, $now[4]+1, $now[3],$now[2],$now[1],$now[0]);
	print STDERR "Timestamp $timestamp\n";
	if (!(-e "$outputDirectory\/metadata")){
	    `mkdir $outputDirectory\/metadata`;
	}
	if (!(-e "$outputDirectory\/$project_name")){
	    `mkdir $outputDirectory\/$project_name`;
	}
	open(VER,">$outputDirectory\/metadata\/Ceto.run.$type.$timestamp.txt");
	&datePrint("Creating metadata file for project $project_name");
	print VER "REF $reference{$project_name}\n";
	print VER "REF Bowtie Index: $bowtieIndex{$reference{$project_name}}\n";
	print VER "REF BWA Index: $bwaIndex{$reference{$project_name}}\n";
	print VER "REF STAR Index: $starIndex{$reference{$project_name}}\n";
	print VER "REF Transcriptome Index: $txIndex{$reference{$project_name}}\n";
	print VER "REF TXDB file: $txdbfile{$reference{$project_name}}\n";
	print VER "REF Exon bed file: $exonbed{$reference{$project_name}}\n";
	print VER "REF GTF file: $gff{$reference{$project_name}}\n";
	print VER "REF Gene bed file: $genebed{$reference{$project_name}}\n";
	print VER "REF RSEM tx file: $rsemTx{$reference{$project_name}}\n";
	print VER "REF GATK reference: $gatkRef{$reference{$project_name}}\n";
	print VER "REF Known SNP sites: $knownSNPsites{$reference{$project_name}}\n";
	print VER "REF Known Indel sites: $knownIndelsites{$reference{$project_name}}\n";
    } 
    if ($fastqDirectory ne ""){
	print VER "INPUT $project_name @fastqlist\n";
	foreach my $fastq (@fastqlist){
	    #	    print STDERR "Fastq: \"$fastq\"\n";
	    if ($fastq ne ""){
		if (($fastq =~ /\/?([\w\-\d\_\.]+)\_S\d/) ||
		    #		($fastq =~ /\/?([\w\-\d\_\.]+)FastqRd/) ||
		    ($fastq =~ /\/?([\w\-\d\_\.]+)\_R\d/) ||
		    ($fastq =~ /\/?([\w\-\d\_\.]+SRR\w+)\_\d.fastq.gz/) ||
		    ($fastq =~ /\/?([\w\-\d\_\.]+).fastq.t?gz/) ||
		    ($fastq =~ /\/?([\w\-\d\_\.]+).fq.t?gz/)
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
	}
    } elsif ($bamDirectory ne ""){
	if (!(-e "$outputDirectory\/$project_name\/scripts")){
	    `mkdir $outputDirectory\/$project_name\/scripts`;
	}
	if (!(-e "$outputDirectory\/$project_name\/bam")){
	    `mkdir $outputDirectory\/$project_name\/bam`;
	}
	&datePrint("Copying bam files into $outputDirectory\/$project_name\/bam\/");
	$cmd .= "cp $bamDirectory\/*.bam $outputDirectory\/$project_name\/bam\/\n";
	system($cmd);
	$bamDirectory = "$outputDirectory\/$project_name\/bam\/";
	$bamlist = `ls $bamDirectory\/*.bam`;
	@bamlist = split(/\s+/,$bamlist);
	print VER "INPUT $project_name @bamlist\n";
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

#print STDERR "Samples list\n";
#print STDERR keys(%samples);

if (($type eq "4C") && ($build4C == 1)){
    my $maxPrimerMismatch = 1;
    if ($fourCdescription ne ""){
	&datePrint("4C experiment described in $fourCdescription");
	$cmd = "mkdir $outputDirectory\/$sample_project\n";
	$cmd .= "mkdir $outputDirectory\/$sample_project\/scripts\n";
	$cmd .= "mkdir $outputDirectory\/$sample_project\/fastq\n";
	system($cmd);
	open(FCSH,">$outputDirectory\/$sample_project\/scripts\/run_4C_demultiplex.sh");
	my (%modulesLoaded, $moduleText);
	print FCSH "$header";
	if ($scheduler eq "MOAB"){
	    print FCSH "#MSUB -N 4Cdemultiplex\n";
	    print FCSH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	} elsif ($scheduler eq "SLURM"){
	    print FCSH "#SBATCH --job-name=4Cdemultiplex\n";
	    print FCSH "#SBATCH --nodes=1\n";
	    print FCSH "#SBATCH -n $numProcessors\n";
	}
	print FCSH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	print FCSH "export PATH=\$PATH:$NGSbartom/tools/bin/\n";
	print FCSH "date\n";
	print FCSH "\n# Copy raw reads into fastq directory and de-compress them.\n";
	print FCSH "cp $baseSpaceDirectory\/Data\/Intensities\/BaseCalls\/Undetermined*.fastq.gz $outputDirectory\/$sample_project\/fastq\/\n";
	print FCSH "gunzip $outputDirectory\/$sample_project\/fastq\/Undetermined*fastq.gz\n";
	print FCSH "date\n";
	print FCSH "\n# De-multiplex reads using $fourCdescription table\n";
	print FCSH "perl $NGSbartom/tools/bin/split4CwithTable.pl $outputDirectory\/$sample_project\/fastq\/ $fourCdescription $maxPrimerMismatch\n";
	print FCSH "date\n";
	print FCSH "\n# Clean up extra files.\n";
	print FCSH "ls $outputDirectory\/$sample_project\/fastq\/*\_S*\_L00*\_R*fastq\n";
	print FCSH "rm $outputDirectory\/$sample_project\/fastq\/*\_S*\_L00*\_R*fastq\n";
	print FCSH "date\n";
	close(FCSH);
	if (($run4C ==1)){
	    my $result = "";
	    &datePrint("Submitting job to de-multiplex 4C samples");
	    if ($scheduler eq "MOAB"){
		$result = `msub $outputDirectory/$sample_project/scripts/run_4C_demultiplex.sh`;
		$result =~ s/\s//g;
	    } elsif ($scheduler eq "SLURM"){
		$result = `sbatch $outputDirectory/$sample_project/scripts/run_4C_demultiplex.sh`;
		$result =~ s/Submitted batch job //g;
		$result =~ s/\s//g;
#		print STDERR "Predicted SLURM job ID:  $result\n";
	    }
	    my $jobfinished = "no";
	    # Wait until the job is Complete.
	    &datePrint("Waiting for job $result to finish. (each . = 300 seconds)");
	    # Check qstat every 300 seconds, adding a "." every time you check.
		waitForJob($result);
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
&datePrint("Setting up directory structure and maybe moving fastq or bam files to project sub-directory within $outputDirectory");
foreach my $project (keys(%samples)){
    if ($scientists{$project} ne ""){
	$scientist = $scientists{$project};
#	if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	$s3path = "ash-tracks/TANGO/$scientist";
    }
    my $cmd = "mkdir $outputDirectory\/$project\n";
    $cmd .= "mkdir $outputDirectory\/$project/scripts\n";
    $cmd .= "mkdir $outputDirectory\/$project\/bam\n";
    $cmd .= "mkdir $outputDirectory\/$project\/fastq\n";
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

# If the aligner is tophat or star, create the shell scripts to run the aligner on all fastqs, one sample at a time.
if (($buildAlign == 1) && ($type eq "RNA")){
    &datePrint("Creating $aligner Alignment shell scripts for RNA-seq analysis.");
    my @samples;
    # Foreach project (TANGO/MOLNG):
    foreach my $project (keys(%samples)){
	if ($scientists{$project} ne ""){
	    $scientist = $scientists{$project};
	    #if ($s3path eq "ash-tracks/TANGO/XXX"){ }
	    $s3path = "ash-tracks/TANGO/$scientist";
	}
	if ($aligner eq "tophat"){
	    # Make a directory for the output.
	    $cmd = "mkdir $outputDirectory\/$project\/Tophat_aln";
	    system($cmd);
	}
	if ($aligner eq "star"){
	    # Make a directory for the output.
	    $cmd = "mkdir $outputDirectory\/$project\/STAR_aln";
	    system($cmd);
	}
	@samples = uniq(split(/\,/,$samples{$project}));

	# Foreach sample within the project:
	foreach my $sample (@samples){
	    # Create a shell script to run tophat on all fastqs for the sample at the same time.
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_$aligner\_align.sh";
	    &datePrint("Printing to $shScript");
	    open (SH,">$shScript");
	    my (%modulesLoaded, $moduleText);
	    print SH "$header";
	    if ($scheduler eq "MOAB"){
		print SH "#MSUB -N $sample\_$aligner\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    } elsif ($scheduler eq "SLURM"){
		print SH "#SBATCH --job-name=$sample\_$aligner\n";
		print SH "#SBATCH --nodes=1\n";
		print SH "#SBATCH -n $numProcessors\n";
	    }
	    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	    print SH "export PATH=\$PATH:$NGSbartom/tools/bin/\n";
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	    #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	    #if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "module load boost/1.56.0\n";
	    print VER "EXEC module load boost/1.56.0\n";
	    print SH "module load gcc/4.8.3\n";
	    print VER "EXEC module load gcc/4.8.3\n";

	    my $rgString = "";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    my $cmd = "";
	    if (!(-e "$outputDirectory\/$project\/fastq")){
		$cmd .= "mkdir $outputDirectory\/$project\/fastq\n";
	    }
	    if (!(-e "$outputDirectory\/$project\/fastqc")){
		$cmd .= "mkdir $outputDirectory\/$project\/fastqc\n";
	    }
	    my $oldfastqs = "@fastqs";
	    my @newfastqs = "";
	    for my $fastq (@fastqs){
		my $basefilename = basename($fastq);
		$fastq = "$outputDirectory\/$project\/fastq\/$basefilename";
		push (@newfastqs, $fastq);
	    }
	    my $newfastqs = "@newfastqs";
	    #	    $newfastqs =~ s/ /,/g;
	    if ($oldfastqs ne $newfastqs){
		print STDERR "Old fastqs: $oldfastqs\nNew fastqs: $newfastqs\n";
		$cmd .= "cp $oldfastqs $outputDirectory\/$project\/fastq/\n";
		system($cmd);
	    }
	    @fastqs = split(/\s+/,$newfastqs);
	    $newfastqs =~ s/\s+/\,/g;
	    if ($newfastqs =~ /^\,([\w\-\_\,\/\.]+)$/){ $newfastqs = $1;}
	    $fastqs{$sample} = $newfastqs;
	    if ($runPairedEnd == 0){
		if ($runTrim == 1){
		    &datePrint("Setting up trimming for single end reads.");
		    print SH "\n# Setting up trimming for single end reads.\n";
# Pascal says it is not necessary to explicitly load java (2018-07-05)
#		    print SH "module load java/jdk1.8.0_25\n";
#		    print VER "EXEC module load java/jdk1.8.0_25\n";
		    my @newfastqs = ();
		    foreach my $fastq (@fastqs){
			if ($fastq ne ""){
			    my $fastqname = "";
			    my $newfastq = "";
			    if (($fastq =~ /\/?([\w\d\-\_\.]+\.fastq\.t?gz$)/) || ($fastq =~ /\/?([\w\d\-\_\.]+\.fq\.t?gz$)/)){
				$fastqname = $1;
			    } elsif ( ($fastq =~ /\/?([\w\d\-\_\.]+\.fastq$)/) || ($fastq =~ /\/?([\w\d\-\_\.]+\.fastq$)/)){
				$fastqname = $1;
			    }
			    $newfastq = "$outputDirectory\/$project\/fastq\/$fastqname";
			    if (!(-e $newfastq) || (-z $newfastq)){
				&datePrint("Copying $fastq to $newfastq\n");
				`cp $fastq $newfastq`;
				$fastq = $newfastq;
			    }
			    print STDERR "Fastq: $fastq\nNewFastq: $newfastq\nFastqname = $fastqname\n";
			    print SH "\n# Trim SE poor quality sequence with $trimString (see Trimmomatic documentation)\n";
			    print SH "java -jar $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $outputDirectory\/$project\/fastq\/$fastqname.trimmed $trimString\n";
			    print VER "EXEC $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar\n";
			    print SH "gzip $outputDirectory\/$project\/fastq\/$fastqname.trimmed\n";
			    if ($newfastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.t?gz$/){
				if ((-e "$1\_fastqc.html") && !(-z "$1\_fastqc.html")){
				    print SH "\# FastQC file already exists\n";
				} else {
				    print SH "# Running FastQC to assess read quality.\n";
				    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $fastq $fastq.trimmed.gz\n";
				    print VER "EXEC $NGSbartom/tools/bin/FastQC/fastqc (FastQC v0.11.2)\n";
				}
				print SH "\n# Running FastQ_screen to look for contamination.\n";
				$moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
				if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
				$moduleText = &checkLoad("perl/5.16",\%modulesLoaded);
				if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"perl/5.16"} = 1;}
				print SH "# Check all reference genomes currently installed\n";
				print SH "perl $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen --threads $numProcessors --aligner bowtie2 --conf $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir $outputDirectory\/$project\/fastqc $fastq\n";
				print VER "EXEC $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen\n";
			    }
			    #		    print SH "mv $newfastq.trimmed.gz $newfastq\n";
			    print SH "date\n\n";
			    push(@newfastqs,"$newfastq.trimmed.gz");
			}
		    }
		    $fastqs{$sample} = "@newfastqs";
		    $fastqs{$sample} =~ s/\s/\,/g;
		    #		&datePrint("New fastqs for sample $sample are @newfastqs, and $fastqs{$sample}");
		}
		if ($genomeBAM == 1){
		    if ($aligner eq "tophat"){
			print SH "\n# Run Tophat to align data for single end data.\n";
			# Add read groups.  These could be made more accurate, by linking into actual metadata rather than an estimate.
			if ($rgString eq ""){
			    $rgString = "--rg-sample $sample --rg-id $sample --rg-library $sample --rg-description $sample --rg-platform-unit nextseq --rg-center ASH --rg-platform nextseq";
			}
			print SH "# Adding Readgroups from rgstring $rgString\n";
			my $libraryType = "";
			if ($stranded == 0) { $libraryType = "--library-type fr-unstranded";}
			$moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
			$moduleText = &checkLoad("tophat/2.1.0",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"tophat/2.1.0"} = 1;}
			if ($reference{$sample} =~ /\.mp/){
			    print SH "\ntophat -G $gff{$reference{$sample}} --read-mismatches $tophatReadMismatch --read-edit-dist $tophatReadEditDist --num-threads $numProcessors $rgString $libraryType -o $outputDirectory\/$project\/Tophat_aln\/$sample $bowtieIndex{$reference{$sample}} $fastqs{$sample} >& $outputDirectory\/$project\/bam\/$sample.tophat.log\n";
			}else {
			    print SH "\ntophat --no-novel-juncs --read-mismatches $tophatReadMismatch --read-edit-dist $tophatReadEditDist --num-threads $numProcessors --max-multihits $tophatMultimap $rgString $libraryType --transcriptome-index $txIndex{$reference{$sample}} -o $outputDirectory\/$project\/Tophat_aln\/$sample $bowtieIndex{$reference{$sample}} $fastqs{$sample} >& $outputDirectory\/$project\/bam\/$sample.tophat.log\n";
			}
		    } elsif ($aligner eq "star"){
			print SH "# Align reads with STAR\n";
			$moduleText = &checkLoad("STAR/2.5.2",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"STAR/2.5.2"} = 1;}
			if ($rgString eq ""){
			    my $date = $timestamp2;
			    if ($sample =~ /\_(\d+)/){
				$date = $1;
			    }
			    $rgString = "ID:$sample LB:$sample PU:nextseq DT:$date SM:$sample CN:ASH PL:illumina";
			}
			print SH "# Adding Readgroups from rgstring $rgString\n";
			print SH "STAR --runMode alignReads --genomeDir $starIndex{$reference{$sample}} --runThreadN $numProcessors --readFilesIn $fastqs{$sample} --readFilesCommand zcat -c --outFileNamePrefix $outputDirectory\/$project\/STAR_aln\/$sample --outSAMtype BAM Unsorted --chimSegmentMin 20 --quantMode TranscriptomeSAM --outReadsUnmapped Fastq --outMultimapperOrder Random --outSAMmapqUnique 60 --outSAMattrRGline $rgString --outFilterMultimapNmax $tophatMultimap --outFilterMismatchNmax $tophatReadMismatch\n\n";
			# Load Samtools if it is not already loaded.
			$moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
			#$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
			#if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
			print SH "# Sort the output of STAR (outputting sorted BAMs from STAR took too much memory)\n";
			print SH "samtools sort -o $outputDirectory\/$project\/STAR_aln\/$sample"."Aligned.sortedByCoord.out.bam $outputDirectory\/$project\/STAR_aln\/$sample"."Aligned.out.bam\n";
		    }
		}
	    } elsif ($runPairedEnd == 1){
		my @read1fastqs;
		my @read2fastqs;
		my @read3fastqs;
		my @fastqs = split(/\,/,$fastqs{$sample});
		&datePrint("All Fastqs for Sample $sample: \n@fastqs");
		if ($fastqs{$sample} !~ /$outputDirectory\/$project\/fastq/){
		    &datePrint("Copying Fastqs to $outputDirectory\/$project\/fastq");
		    `cp @fastqs $outputDirectory\/$project\/fastq/`;
		}
		foreach my $fastq (@fastqs){
		    #		    print STDERR "Looking for Read number in $fastq\n";
		    if (($fastq =~ /\_R1\_?\.?/) || ($fastq =~ /SRR.+\_1\./) ){
			#			|| ($fastq =~ /Rd1/)){
			push (@read1fastqs,$fastq);
		    } elsif (($fastq =~ /\_R2\_?\.?/) || ($fastq =~ /SRR.+\_2\./)){
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
		print STDERR "Read1fastqs: $read1fastqs\n";
		print STDERR "Read2fastqs: $read2fastqs\n";
		print STDERR "Read3fastqs: $read3fastqs\n";
		if (length($read3fastqs)>length($read2fastqs)){
		    $read2fastqs = $read3fastqs;
		}
		$read1fastqs =~ s/\ /,/g;
		$read2fastqs =~ s/\ /,/g;
		print SH "\n# Running FastQ_screen to look for contamination.\n";
		$moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
		$moduleText = &checkLoad("perl/5.16",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"perl/5.16"} = 1;}
		print SH "# Check all reference genomes currently installed\n";
		print SH "perl $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen --threads $numProcessors --aligner bowtie2 --conf $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir $outputDirectory\/$project\/fastqc $read1fastqs\n";
		print VER "EXEC $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen\n";

		if ($runTrim == 1){
		    &datePrint("Setting up trimming for sample $sample paired end reads.");
		    print SH "\n# Setting up trimming for sample $sample for paired end reads.\n";
		    print SH "\n# Trim PE poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads $numProcessors -phred33 $read1fastqs $read2fastqs $outputDirectory/$project/fastq/$sample\_R1.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R1U.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R2.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R2U.fastq.trimmed.gz $trimString\n\n";
		    print VER "EXEC $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar\n";
		    print SH "# Running FastQC to assess read quality.\n";
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $read1fastqs $read2fastqs\n";
		    print VER "EXEC $NGSbartom/tools/bin/FastQC/fastqc (FastQC v0.11.2)\n";
		    my $trimmedUnpaired1 = $read1fastqs;
		    my $trimmedUnpaired2 = $read2fastqs;
		    $trimmedUnpaired1 =~ s/R1.fastq/R1U.fastq.trimmed/g;
		    $trimmedUnpaired2 =~ s/R2.fastq/R2U.fastq.trimmed/g;
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $trimmedUnpaired1 $trimmedUnpaired2\n";
		    my $trimmedPaired1 = $read1fastqs;
		    my $trimmedPaired2 = $read2fastqs;
		    $trimmedPaired1 =~ s/R1.fastq/R1.fastq.trimmed/g;
		    $trimmedPaired2 =~ s/R2.fastq/R2.fastq.trimmed/g;
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $trimmedPaired1 $trimmedPaired2\n";
		    print SH "date\n\n";
		    $read1fastqs = $trimmedPaired1;
		    $read2fastqs = $trimmedPaired2;
		}
		if ($genomeBAM == 1){
		    if ($aligner eq "tophat"){
			print SH "\n# Run Tophat to align data for paired end data.\n";
			my $libraryType = "";
			if ($stranded == 0) { $libraryType = "--library-type fr-unstranded";}
			if ($rgString eq ""){
			    $rgString = "--rg-sample $sample --rg-id $sample --rg-library $sample --rg-description $sample --rg-platform-unit nextseq --rg-center ASH --rg-platform nextseq";
			}
			$moduleText = &checkLoad("tophat/2.1.0",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"tophat/2.1.0"} = 1;}
			print SH "# Adding Readgroups from rgstring $rgString\n";
			print SH "\ntophat --no-novel-juncs --read-mismatches $tophatReadMismatch --read-edit-dist $tophatReadEditDist --num-threads $numProcessors --max-multihits $tophatMultimap $rgString $libraryType --transcriptome-index $txIndex{$reference{$sample}} -o $outputDirectory\/$project\/Tophat_aln\/$sample $bowtieIndex{$reference{$sample}} $read1fastqs $read2fastqs >& $outputDirectory\/$project\/bam\/$sample.tophat.log\n";
		    } elsif ($aligner eq "star"){
			print SH "# Align reads with STAR\n";
			$moduleText = &checkLoad("STAR/2.5.2",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"STAR/2.5.2"} = 1;}
			if ($rgString eq ""){
			    my $date = $timestamp2;
			    if ($sample =~ /\_(\d+)/){
				$date = $1;
			    }
			    $rgString = "ID:$sample LB:$sample PU:nextseq DT:$date SM:$sample CN:ASH PL:illumina";
			}
			print SH "# Adding Readgroups from rgstring $rgString\n";
			print SH "STAR --runMode alignReads --genomeDir $starIndex{$reference{$sample}} --runThreadN $numProcessors --readFilesIn $read1fastqs $read2fastqs --readFilesCommand zcat -c --outFileNamePrefix $outputDirectory\/$project\/STAR_aln\/$sample --outSAMtype BAM Unsorted --chimSegmentMin 20 --quantMode TranscriptomeSAM --outReadsUnmapped Fastq --outMultimapperOrder Random --outSAMattrRGline $rgString --outSAMmapqUnique 60 --outFilterMultimapNmax $tophatMultimap --outFilterMismatchNmax $tophatReadMismatch\n\n";
			print SH "# Sort the output of STAR (outputting sorted BAMs from STAR took too much memory)\n";
			$moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
			if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
			#$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
			#if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
			print SH "samtools sort -o $outputDirectory\/$project\/STAR_aln\/$sample"."Aligned.sortedByCoord.out.bam $outputDirectory\/$project\/STAR_aln\/$sample"."Aligned.out.bam\n";
		    }
		}
		if ($rsem == 1){ # AND runpaired = 1
		    $moduleText = &checkLoad("STAR/2.5.2",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"STAR/2.5.2"} = 1;}
		    $moduleText = &checkLoad("perl/5.16",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"perl/5.1.6"} = 1;}
		    print SH "module load gcc/6.4.0\n";
		    print VER "EXEC module load gcc/6.4.0\n";
		    print SH "export PATH=$NGSbartom/tools/bin/RSEM-1.3.0/:\$PATH\n";
		    print VER "EXEC $NGSbartom/tools/bin/RSEM-1.3.0\n";
		    print SH "date\n";
		    my $strandstring = "";
		    if ($stranded == 1) { $strandstring = "--forward-prob 0";}
		    if ($aligner ne "star"){
			print SH "\n# NOTE: This is buggy. I recommend using the STAR aligner with RSEM.\n";
			print SH "\n# Align fastqs to transcriptome for RSEM\n";
			print SH "\n# First prepare fastqs\n";
			#		print SH "mkdir $outputDirectory/$project/fastq/\n";
			$read1fastqs =~ s/\,/\ /g;
			$read2fastqs =~ s/\,/\ /g;
			&datePrint("Fastqs for RSEM are $read1fastqs and $read2fastqs.");
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
			print SH "which rsem-calculate-expression\n";
			print SH "rsem-calculate-expression --star --paired-end $outputDirectory/$project/fastq/$sample.read1.fastq $outputDirectory/$project/fastq/$sample.read2.fastq $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample --no-bam-output -p $numProcessors $strandstring >& $outputDirectory/$project/bam/$sample.rsem.log\n";
			print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n";
			print SH "ls $outputDirectory/$project/fastq/$sample.read*.fastq\n";
			print SH "rm $outputDirectory/$project/fastq/$sample.read1.fastq\n";
			print SH "rm $outputDirectory/$project/fastq/$sample.read2.fastq\n";
		    }elsif ($aligner eq "star"){
			print SH "\n# Use STAR alignment to transcriptome and RSEM to estimate isoform expression levels.\n";
			print SH "rsem-calculate-expression --alignments --paired-end $outputDirectory/$project/STAR_aln/$sample"."Aligned.toTranscriptome.out.bam $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample -p $numProcessors $strandstring --no-bam-output >& $outputDirectory/$project/bam/$sample.rsem.log\n";
			print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n\n";
		    }
		}
	    }
	    print SH "date\n";
	#	    print SH "\nrsync -av \$TMPDIR\/$sample/* $outputDirectory\/$project\/Tophat_aln\/$sample/\n";
	    #	    print SH "\nmv \$TMPDIR\/$project\_$sample\_tophatOut $outputDirectory\/$project\/Tophat_aln\/$sample\n";
	    print SH "date\n";
	    if ($aligner eq "tophat"){
		print SH "ln -s $outputDirectory\/$project\/Tophat_aln\/$sample\/accepted_hits.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    }
	    if ($aligner eq "star"){
		print SH "ln -s $outputDirectory\/$project\/STAR_aln\/$sample"."Aligned.sortedByCoord.out.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    }
	    print SH "date\n";
	    # Modifications made so that htseq will work correctly on paired end data.
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	    #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	    #if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    if (($runPairedEnd == 1) && ($htseq == 1)) {
		print SH "cd $outputDirectory\/$project\/bam\n";
		print SH "samtools sort -n -T $sample -o $outputDirectory\/$project\/bam\/$sample.sorted.names.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
		print SH "mv $outputDirectory\/$project\/bam\/$sample.sorted.names.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
		print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    } else {
		print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    }
	    if (($rsem == 1)&& ($runPairedEnd == 0)){ # Single End RSEM analysis
		print SH "\n# Use RSEM to analyze isoform abundance.\n";
		print SH "export PATH=$NGSbartom/tools/bin/RSEM-1.3.0/:\$PATH\n";
		print VER "EXEC module load gcc/6.4.0\n";
		print SH "module load gcc/6.4.0\n";
		$moduleText = &checkLoad("perl/5.16",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"perl/5.16"} = 1;}
		$moduleText = &checkLoad("STAR/2.5.2",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"STAR/2.5.2"} = 1;}
		print VER "EXEC $NGSbartom/tools/bin/RSEM-1.3.0\n";
		print SH "date\n";
		my $strandstring = "";
		if ($stranded == 1) { $strandstring = "--forward-prob 0";}
		if ($aligner ne "star"){
		    print SH "\n# NOTE: This is buggy. I recommend using the STAR aligner with RSEM.\n";
		    print SH "\n# First align fastqs to transcriptome with STAR\n";
		    my $fastqstring = $fastqs{$sample};
		    $fastqstring =~ s/\,/ /g;
		    print SH "gunzip $fastqstring\n";
		    $fastqstring =~ s/.gz//g;
		    print SH "echo $fastqstring\n";
		    print SH "cat $fastqstring > $outputDirectory/$project/fastq/$sample.fastq\n";
		    print SH "ls $outputDirectory/$project/fastq/$sample.fastq\n";
		    print SH "gzip $fastqstring &\n";
		    print SH "\n# Then calculate expression for $sample.\n";
		    $moduleText = &checkLoad("STAR/2.5.2",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"STAR/2.5.2"} = 1;}
		    print SH "rsem-calculate-expression --star $outputDirectory/$project/fastq/$sample.fastq $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample --no-bam-output -p $numProcessors $strandstring >& $outputDirectory/$project/bam/$sample.rsem.log\n";
		    print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n";
		    print SH "ls $outputDirectory/$project/fastq/$sample.fastq\n";
		    print SH "# rm $outputDirectory/$project/fastq/$sample.fastq\n";
		}elsif ($aligner eq "star"){
		    print SH "\n# Use STAR alignment to transcriptome and RSEM to estimate isoform expression levels.\n";
			print SH "rsem-calculate-expression --alignments $outputDirectory/$project/STAR_aln/$sample"."Aligned.toTranscriptome.out.bam $rsemTx{$reference{$sample}} $outputDirectory/$project/bam/$sample --no-bam-output -p $numProcessors $strandstring  >& $outputDirectory/$project/bam/$sample.rsem.log\n";
		    print SH "rsem-plot-model $outputDirectory/$project/bam/$sample $outputDirectory/$project/bam/$sample.rsemPlot.pdf\n\n";
		}
	    }
	    if ($htseq == 1) {
		print SH "module unload mpi\n";
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
		print VER "EXEC htseq 0.6.1\n";
		print SH "\n# Run htseq-count for $sample, stranded = $stranded\n";
		if ($stranded == 1){
		    print SH "htseq-count -f bam -q -m intersection-nonempty -s reverse -t exon -i gene_id $outputDirectory\/$project\/bam\/$sample.bam $gff{$reference{$sample}} > $outputDirectory\/$project\/bam\/$sample.htseq.counts\n";
		} elsif ($stranded == 0){
		    print SH "htseq-count -f bam -q -m intersection-nonempty -s no -t exon -i gene_id $outputDirectory\/$project\/bam\/$sample.bam $gff{$reference{$sample}} > $outputDirectory\/$project\/bam\/$sample.htseq.counts\n";
		}
		if ($modulesLoaded{"python/anaconda"} == 1){ print SH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
		print SH "module load gcc/4.8.3\n";
	    }
	    if ($bedtools == 1) {
		$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
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
		$modulesLoaded{"bedtools/2.17.0"} = 0;
	    }
	    print SH "date\n\n";
	    if ($makeTracks == 1){
		print SH "# Create RNA seq Tracks\n";
		print SH "module load R/3.2.2\n";
		print VER "EXEC module load R/3.2.2\n";
		if ($stranded == 1){
		    my $assembly = $reference{$sample};
		    if ($assembly eq "hg38.mp"){$assembly="hg38";}
		    print SH "Rscript $NGSbartom/tools/bin/createRNAseqTracks3.R --assembly=$assembly --bamDir=$outputDirectory\/$project\/bam\/ --sample=$sample --multiMap=$multiMap --stranded=$stranded\n";
		} else {
		    my $assembly = $reference{$sample};
		    if ($assembly eq "hg38.mp"){$assembly="hg38";}
		    # The multi mapping argument is used here, because these are mapped with tophat.
		    #		    print SH "Rscript $NGSbartom/tools/bin/createChIPtracks2.R --assembly=$assembly --bamDir=$outputDirectory\/$project\/bam\/ --sample=$sample --extLen=0 --multiMap=$multiMap\n";
		    print SH "Rscript $NGSbartom/tools/bin/createRNAseqTracks3.R --assembly=$assembly --bamDir=$outputDirectory\/$project\/bam\/ --sample=$sample --multiMap=$multiMap --stranded=$stranded\n";
		}
	    print SH "date\n\n";
		print SH "mkdir $outputDirectory\/$project\/tracks\n";
		print SH "\n# Make Headers for UCSC genome browser.\n";
		if ($stranded == 1){
		    if ($multiMap == 0) {
			print SH "echo \"track type=bigWig name=$sample.plus.bw description=$sample.plus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=255,0,0 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.plus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.plus.bw.header.txt\n";
			print SH "echo \"track type=bigWig name=$sample.minus.bw description=$sample.minus.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.minus.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.minus.bw.header.txt\n";
		    } elsif ($multiMap == 1){
			print SH "echo \"track type=bigWig name=$sample.plus.multi.bw description=$sample.plus.multi.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=255,0,0 itemRGB=on autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.plus.multi.bw\" | cat > $outputDirectory\/$project\/tracks\/$sample.plus.bw.header.multi.txt\n";
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
#		print SH "\nmodule load python/anaconda\n";

		print SH "# Move tracks into Shilatifard directory structure.\n";
		print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
		    print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
		print SH "cp $outputDirectory\/$project\/bam\/$sample*.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
		    if ($uploadBAM == 1){
			print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "# To load, paste the following url into the custom tracks field:\n";
			print SH "# http://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/$scientist/$scientist.$project/$sample.bam\n";
			if (($runPairedEnd == 1) && ($htseq == 1)){
			    print SH "# If files are paired end and htseq is used for gene counting, then BAMs are sorted by name and indices are not created.  This will require an extra re-sorting and indexing step, to be added if this is a common concern.\n";
			}
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://$s3path/$scientist.$project/ --region us-west-2\n";
			print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://$s3path/$scientist.$project/ --region us-west-2\n";
		    }
		    print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
		    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
		    if ($stranded == 1){
			if ($multiMap == 0){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
		    } else {
			if ($multiMap == 0){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}elsif ($multiMap == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.multi.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
		    }
		if ($modulesLoaded{"python/anaconda"} == 1){ print SH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
	    } elsif ($uploadPulmtracks == 1) {
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
		if ($uploadBAM == 1){
		    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
		    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
		    print SH "# To load, paste the following url into the custom tracks field:\n";
		    print SH "# http://s3-us-west-2.amazonaws.com/m-328-data/$scientist/$scientist.$project/$sample.bam\n";
		    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
		    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
		}
	    } elsif ($uploadLutracks == 1) {
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
		if ($uploadBAM == 1){
		    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
		    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
		    print SH "# To load, paste the following url into the custom tracks field:\n";
		    print SH "# http://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$sample.bam\n";
		    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://lwwlab/$scientist.$project/ --region us-west-2\n";
		    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://lwwlab/$scientist.$project/ --region us-west-2\n";
		}
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
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
		if ($modulesLoaded{"python/anaconda"} == 1){ print SH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
	    }
	    close(SH);
	}
	if ($runAlign == 0){
	    # Print tips on running the tophat shell scripts.
	    print STDERR "To execute all scripts, use the following command:\n";
	    #	"$outputDirectory\/$project\/scripts\/run\_$sample\_$aligner\_align.sh";
	    if ($scheduler eq "MOAB"){
		print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec msub {} ./ \\\;\n";
	    } elsif ($scheduler eq "SLURM"){
		print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec sbatch {} ./ \\\;\n";
	    }
	}
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
	my $cmd = "mkdir $outputDirectory\/$project\/bam\n";
	$cmd .= "mkdir $outputDirectory\/$project\/fastqc\n";
	system($cmd);
	@samples = uniq(split(/\,/,$samples{$project}));
	print STDERR "Samples: @samples\n";
	# Foreach sample within the project:
	foreach my $sample (@samples){
	    # Create a shell script to run bowtie on all fastqs for the sample at the same time.
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_$aligner\_align.sh";
	    &datePrint("Printing $shScript");
	    open (SH,">$shScript");
	    my (%modulesLoaded, $moduleText);
	    print SH "$header";
	    if ($scheduler eq "MOAB"){
		print SH "#MSUB -N $sample\_$aligner\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    } elsif ($scheduler eq "SLURM"){
		print SH "#SBATCH --job-name=$sample\_$aligner\n";
		print SH "#SBATCH --nodes=1\n";
		print SH "#SBATCH -n $numProcessors\n";
	    }
	    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	    $moduleText = &checkLoad("bowtie/1.1.2",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie/1.1.2"} = 1;}
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	   #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	    #if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "module load R/3.2.2\n";
	    print VER "EXEC module load R/3.2.2\n";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    my @newFastqs;
	    if ($runTrim == 1){
		# Pascal says this is unnecessary as java is loaded by default.  2018-07-05
#		print SH "module load java/jdk1.8.0_25\n";
		foreach my $fastq (@fastqs){
		    print SH "\n# Copy fastq files into outputDirectory.\n";
		    print SH "cp $fastq $outputDirectory/$project/fastq/\n";
		    if ($fastq =~ /\/?([\d\_\-\w\.\.]+)\.fastq\.gz$/){
			$fastq = "$outputDirectory/$project/fastq/$1.fastq.gz";
			push(@newFastqs,$fastq);
		    }
		    print SH "\n# Trim SE poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $fastq.trimmed $trimString\n";
		    print VER "EXEC  $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar\n";
		    print SH "gzip $fastq.trimmed\n";
		    if ($fastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.gz$/){
			if (-f "$1\_fastqc.html"){
			    print SH "\# FastQC file already exists\n";
			} else {
			    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc $fastq $fastq.trimmed.gz\n";
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
	    if ($runTrim == 0){
		foreach my $fastq (@fastqs){
		    print SH "\n# Copy fastq files into outputDirectory.\n";
		    print SH "cp $fastq $outputDirectory/$project/fastq/\n";
		    if ($fastq =~ /\/?([\d\_\-\w\.]+).fastq\.gz$/){
			$prefix = $1;
		    } else { $prefix = $fastq;}
		    $fastq = "$outputDirectory\/$project/fastq/$prefix.fastq.gz";
		    push(@newFastqs,$fastq);
		}
	    }
	    @fastqs = @newFastqs;
	    my $fastqs = "@fastqs";
	    $fastqs =~ s/,/ /g;
	    if ($runTrim == 1){
		$fastqs =~ "s/q.gz/q.trimmed.gz/g";
	    }
	    print SH "date\n\n";
	    #	    print STDERR "$sample\tREF:$reference{$sample}\tINDEX:$bowtieIndex{$reference{$sample}}\n";
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	    #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	    #if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "# Align fastqs with Bowtie\n";
	    if ($multiMap == 0){
		print SH "\ngunzip -c $fastqs | bowtie -p $numProcessors -m 1 -v 2 -S $bowtieIndex{$reference{$sample}} 2> $outputDirectory\/$project\/bam\/$sample.bowtie.log - | samtools view -bS - > $outputDirectory\/$project\/bam\/$sample.bam \n";
	    } elsif ($multiMap ==1) {
		print SH "\ngunzip -c $fastqs | bowtie -p $numProcessors -v 2 -S $bowtieIndex{$reference{$sample}} 2> $outputDirectory\/$project\/bam\/$sample.bowtie.log - | samtools view -bS - > $outputDirectory\/$project\/bam\/$sample.bam \n";
	    }
	    print SH "date\n\n";
	    print SH "# Sort and rearrange bam files\n";
	    print SH "samtools sort -o $outputDirectory\/$project\/bam\/$sample.sorted.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n";
	    print SH "mv $outputDirectory\/$project\/bam\/$sample.sorted.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "date\n";
	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
	    print SH "\n# Check if bamlist file exists, and if not, create it.\n";
	    print SH "if ! [ -f \"$outputDirectory\/$project\/bam\/bamlist.txt\" ];\nthen\necho \"$outputDirectory\/$project\/bam\/$sample.bam\" | cat > $outputDirectory\/$project\/bam\/bamlist.txt \n";
	    print SH "else\necho \"$outputDirectory\/$project\/bam\/$sample.bam\" | cat >> $outputDirectory\/$project\/bam\/bamlist.txt \nfi\n";
	    print SH "date\n\n";
	    if ($type eq "chipseq"){
		if ($makeTracks == 1){
		    print SH "# Make ChIPseq tracks.\n";
		    # The multi mapping argument is not used right now.
		    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
		    my $assembly = $reference{$sample};
		    if ($assembly eq "hg38.mp"){ $assembly = "hg38";}
		    print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$assembly --bamDir=$bamDirectory --sample=$sample --extLen=150\n";
		    print SH "date\n\n";
		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
		    if ($uploadASHtracks == 1){
			print SH "# Move tracks into Shilatifard directory structure\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
			print SH "cp $outputDirectory\/$project\/bam\/$sample*.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			if ($uploadBAM == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
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
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
		    }
		    if ($uploadPulmtracks == 1){
			if ($uploadBAM == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/m-328-data/$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}
			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "aws s3 cp $bamDirectory\/$sample.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
		    }
		    if ($uploadLutracks == 1){
			if ($uploadBAM == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/lwwlab/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://lwwlab/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://lwwlab/$scientist.$project/ --region us-west-2\n";
			}
			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			print SH "aws s3 cp $bamDirectory\/$sample.bw s3://lwwlab/$scientist.$project/ --region us-west-2\n";
		    }
		    print SH "mv $bamDirectory\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
		    print SH "\n# Check if bwlist file exists, and if not, create it.\n";
		    print SH "if ! [ -f \"$outputDirectory\/$project\/tracks\/bwlist.txt\" ];\nthen\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/bwlist.txt \n";
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
			    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
			  #  $moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
			  #  if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
			    print SH "\n# Remove Viewpoint from Bam files.\n";
			    print SH "samtools view -h $bamDirectory\/$sample.bam | awk '!(\$3 == \"$chr\" && \$4 > $start && \$4 < $stop){print \$0}' | samtools view -Sb - > $bamDirectory\/$sample.noVP.bam\n";
			    print SH "\n# Make 4C tracks.\n";
			    # The multi mapping argument is not used right now.
			    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
			    my $assembly = $reference{$sample};
			    if ($assembly eq "hg38.mp"){$assembly = "hg38";}
			    print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$assembly --bamDir=$bamDirectory --sample=$sample.noVP --extLen=0\n";
			    print SH "mv $bamDirectory\/$sample.noVP.bw $bamDirectory\/$sample.bw\n";
			    print SH "date\n";
			}
		    } else {
			print SH "\n# Make 4C tracks.\n";
			# The multi mapping argument is not used right now.
			# This is because Bowtie doesn't fill in the NH tag in the BAM file.
			my $assembly = $reference{$sample};
			if ($assembly eq "hg38.mp"){$assembly = "hg38";}
			print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$assembly --bamDir=$bamDirectory --sample=$sample --extLen=0\n";
		    }
		    print SH "date\n\n";
		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
		    if ($uploadASHtracks == 1){
			print SH "\n# Move tracks into Shilatifard directory structure\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
			print SH "cp $bamDirectory\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
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
	if ($scheduler eq "MOAB"){
	    print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec msub {} ./ \\\;\n";
	} elsif ($scheduler eq "SLURM"){
	    print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec sbatch {} ./ \\\;\n";
	}
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
	    my $shScript = "$outputDirectory\/$project\/scripts\/run\_$sample\_$aligner\_align.sh";
	    &datePrint("Printing $shScript");
	    open (SH,">$shScript");
	    my (%modulesLoaded, $moduleText);
	    print SH "$header";
	    if ($scheduler eq "MOAB"){
		print SH "#MSUB -N $sample\_$aligner\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    } elsif ($scheduler eq "SLURM"){
		print SH "#SBATCH --job-name=$sample\_$aligner\n";
		print SH "#SBATCH --nodes=1\n";
		print SH "#SBATCH -n $numProcessors\n";
	    }
	    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	    $moduleText = &checkLoad("bwa/0.7.12",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bwa/0.7.12"} = 1;}
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	#    $moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	#    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "module load R/3.2.2\n";
	    print VER "EXEC module load R/3.2.2\n";
	    print SH "module load picard/1.131\n";
	    print VER "EXEC module load picard/1.131\n";
	    print SH "module load java/jdk1.8.0_25\n";
	    print VER "EXEC module load java/jdk1.8.0_25\n";
	    print SH "\nmkdir $outputDirectory\/$project\/fastq\n";
	    print SH "\nmkdir $outputDirectory\/$project\/fastqc\n";
	    my $PICARD = "/software/picard/1.131/picard-tools-1.131/picard.jar";
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    my @newFastqs;

	    foreach my $fastq (@fastqs){
		print SH "\n# Running FastQ_screen to look for contamination.\n";
		$moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
		$moduleText = &checkLoad("perl/5.16",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"perl/5.16"} = 1;}
		print SH "# Check all reference genomes currently installed\n";
		print SH "perl $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen --threads $numProcessors --aligner bowtie2 --conf $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen.allRefs.conf --outdir $outputDirectory\/$project\/fastqc $fastq\n";
		print VER "EXEC $NGSbartom/tools/bin/fastq_screen_v0.11.4/fastq_screen\n";
	    }

	    if (($runTrim == 1) && ($runPairedEnd == 0)){
		foreach my $fastq (@fastqs){
		    print SH "\n# Copy single end fastq files into outputDirectory.\n";
		    print SH "cp $fastq $outputDirectory/$project/fastq/\n";
		    if ($fastq =~ /\/?([\d\_\-\w\.\.]+)\.fastq\.gz$/){
			$fastq = "$outputDirectory/$project/fastq/$1.fastq.gz";
			push(@newFastqs,$fastq);
		    }
		    print SH "\n# Trim SE poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar SE -threads $numProcessors -phred33 $fastq $fastq.trimmed $trimString\n";
		    print VER "EXEC  $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar\n";
		    print SH "gzip $fastq.trimmed\n";
		    if ($fastq =~ /^([\d\_\-\w\.\/.]+)\.fastq\.gz$/){
			if (-f "$1\_fastqc.html"){
			    print SH "\# FastQC file already exists\n";
			} else {
			    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc $fastq $fastq.trimmed.gz\n";
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
		    print SH "\n# Align Single End Fastq with BWA\n";
		    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		    $rgString = "\@RG\\\tID:$sample\\\tSM:$sample\\\tPU:nextseq\\\tCN:NUSeq\\\tPL:ILLUMINA";
		    print SH "# Adding Readgroups from rgstring $rgString\n";
		    print SH "bwa mem -M -t $numProcessors -R \"$rgString\" $bwaIndex{$reference{$sample}} $fastq |  $NGSbartom/tools/bin/samblaster/samblaster | samtools view -bS - > $outputDirectory\/$project\/bam\/$prefix.bam\n";
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
		@read1fastqs = uniq(sort(@read1fastqs));
		@read2fastqs = uniq(sort(@read2fastqs));
		@read3fastqs = uniq(sort(@read3fastqs));
		my $read1fastqs = "@read1fastqs";
		my $read2fastqs = "@read2fastqs";
		my $read3fastqs = "@read3fastqs";
		if (length($read3fastqs)>length($read2fastqs)){
		    @read2fastqs = @read3fastqs;
		}
		$read1fastqs =~ s/\ /,/g;
		$read2fastqs =~ s/\ /,/g;
#		print SH "\n# R1: $read1fastqs\n";
		#		print SH "\n# R2: $read2fastqs\n";
		if ($runTrim == 1){
		    &datePrint("Setting up trimming for sample $sample for paired end reads.");
		    print SH "# Setting up trimming for sample $sample for paired end reads.\n";
		    print SH "\n# Trim PE poor quality sequence with $trimString (see Trimmomatic documentation)\n";
		    print SH "java -jar $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads $numProcessors -phred33 $read1fastqs $read2fastqs $outputDirectory/$project/fastq/$sample\_R1.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R1U.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R2.fastq.trimmed.gz $outputDirectory/$project/fastq/$sample\_R2U.fastq.trimmed.gz $trimString\n\n";
		    print VER "EXEC $NGSbartom/tools/bin/Trimmomatic-0.33/trimmomatic-0.33.jar\n";
		    print SH "# Running FastQC to assess read quality before trimming.\n";
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $read1fastqs $read2fastqs\n";
		    print VER "EXEC $NGSbartom/tools/bin/FastQC/fastqc (FastQC v0.11.2)\n";
		    my $trimmedUnpaired1 = "$outputDirectory/$project/fastq/$sample\_R1U.fastq.trimmed.gz";
		    my $trimmedUnpaired2 = "$outputDirectory/$project/fastq/$sample\_R2U.fastq.trimmed.gz";
		    print SH "\n# Running FastQC to assess read quality of unpaired reads after trimming.\n";
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $trimmedUnpaired1 $trimmedUnpaired2\n";
		    my $trimmedPaired1 = "$outputDirectory/$project/fastq/$sample\_R1.fastq.trimmed.gz";
		    my $trimmedPaired2 = "$outputDirectory/$project/fastq/$sample\_R2.fastq.trimmed.gz";
		    print SH "\n# Running FastQC to assess read quality of paired reads after trimming.\n";
		    print SH "date\n$NGSbartom/tools/bin/FastQC/fastqc -o $outputDirectory\/$project\/fastqc $trimmedPaired1 $trimmedPaired2\n";
		    print SH "date\n\n";
		    $read1fastqs = $trimmedPaired1;
		    $read2fastqs = $trimmedPaired2;
		}
		@read1fastqs = split(/\,/,$read1fastqs);
		@read2fastqs = split(/\,/,$read2fastqs);
		# This is somewhat half set up to deal with multiple fastqs for the same sample.  It works if there is only one read1 fastq per sample.  It
		# would need to be updated for two or more read1 fastqs per sample.  #ebartom 2019-09-03
		for (my $i=0;$i <= $#read1fastqs;$i++){
		    my $prefix = "";
		    if ($read1fastqs[$i] =~ /\/?([\d\_\-\w\.]+)\_R1\.fastq\.gz$/){
			$prefix = $1;
		    } elsif ($read1fastqs[$i] =~ /\/?([\d\_\-\w\.]+)\_R1\_001\.fastq\.gz$/){
			$prefix = $1;
		    } else { $prefix = $read1fastqs[$i];}
		    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		    #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
		    #if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
		    print SH "\n# Align Paired End Fastq with BWA\n";
		    #		    print SH "# Line L1909\n";
		    $rgString = "\@RG\\tID:$sample\\tSM:$sample\\tPU:nextseq\\tCN:NUSeq\\tPL:ILLUMINA";
		    print SH "# Adding Readgroups from rgstring $rgString\n";
		    print SH "bwa mem -M -t $numProcessors -R \"$rgString\" $bwaIndex{$reference{$sample}} $read1fastqs[$i] $read2fastqs[$i] | \\\n";
		    print SH "\t$NGSbartom/tools/bin/samblaster/samblaster | \\\n";
		    print SH "\tsamtools view -b -h - | \\\n";
		    print SH "\tsamtools sort -o $outputDirectory\/$project\/bam/$sample.bam - && \\\n";
		    print SH "\tsamtools index $outputDirectory\/$project\/bam\/$sample.bam && \\\n";
		    print SH "\tsamtools flagstat $outputDirectory\/$project\/bam\/$sample.bam && \n";
		    print SH "date\n\n";
		    push (@bams,"$outputDirectory\/$project\/bam\/$sample.bam");
		}
	    }
	    # I wrote this earlier to except multiple read1fastqs for the same BAM, but I don't have an example of that now, so I will assume there is only one
	    # pair of reads per bam; this may need to be fixed later. #ebartom 2019-09-03
#	    my $bamString = "I=@bams";
#	    $bamString =~ s/\s/ I=/g;
#	    print SH "# Merge bam files from the same sample $sample\n";
#	    print SH "java -Xmx2g -jar $PICARD MergeSamFiles $bamString O=$outputDirectory\/$project\/bam\/$sample.bam\n";
	    #	    print SH "date\n\n";
	    print SH "echo \"Finished alignment\"\n";
#	    print SH "# Sort bam file for sample $sample\n";
#	    print SH "java -Xmx2g -jar $PICARD SortSam \\\n";
#	    print SH "\tINPUT=$outputDirectory\/$project\/bam\/$sample.bam \\\n";
#	    print SH "\tOUTPUT=$outputDirectory\/$project\/bam\/$sample.sorted.bam \\\n";
#	    print SH "\tSORT_ORDER=coordinate\n";
	    print SH "date\n\n";
	    print SH "# Mark duplicate reads in the bam file for sample $sample\n";
	    print SH "java -Xmx2g -jar $PICARD MarkDuplicates \\\n";
	    print SH "\tINPUT=$outputDirectory\/$project\/bam\/$sample.bam \\\n";
	    print SH "\tOUTPUT=$outputDirectory\/$project\/bam\/$sample.mdup.bam \\\n";
	    print SH "\tMETRICS_FILE=$outputDirectory\/$project\/bam\/$sample.mdup_metrics.txt\n";
 	    print SH "date\n\n";
	    print SH "# Index bam file and gather flagstats.\n";
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
 	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.mdup.bam\n";
	    print SH "samtools flagstat $outputDirectory\/$project\/bam\/$sample.mdup.bam > $outputDirectory\/$project\/bam\/$sample.mdup.bam.flagstats.txt &\n";
 	    print SH "date\n\n";
	    #Base recalibrator
	    print SH "# Base recalibrator\n";
	    print SH "java -Xmx2G -jar /projects/p20742/tools/bin/GATK_v3.6//GenomeAnalysisTK.jar -T BaseRecalibrator \\\n";
	    print SH "\t-R $bwaIndex{$reference{$project}} \\\n";
	    print SH "\t-I $outputDirectory/$project/bam/$sample.mdup.bam \\\n";
	    print SH "\t-knownSites $knownSNPsites{$reference{$sample}} \\\n"; 
	    print SH "\t-o $outputDirectory/$project/bam/$sample.mdup.recalibration_report.grp\n\n";
	    
	    print SH "# Print Recalibrated base quality scores to bam file\n";
	    print SH "java -Xmx2G -jar /projects/p20742/tools/bin/GATK_v3.6/GenomeAnalysisTK.jar -T PrintReads \\\n";
	    print SH "\t-R $bwaIndex{$reference{$project}} \\\n";
	    print SH "\t-I $outputDirectory/$project/bam/$sample.mdup.bam \\\n";
	    print SH "\t-BQSR $outputDirectory/$project/bam/$sample.mdup.recalibration_report.grp \\\n";
	    print SH "\t-o $outputDirectory/$project/bam/$sample.recal.bam\n\n";

	    print SH "# Calculate depth of coverage\n";
	    print SH "java -Xmx2G -jar /projects/p20742/tools/bin/GATK_v3.6/GenomeAnalysisTK.jar -T DepthOfCoverage --omitDepthOutputAtEachBase \\\n";
	    print SH "\t--summaryCoverageThreshold 10 --summaryCoverageThreshold 25 --summaryCoverageThreshold 50 --summaryCoverageThreshold 100 \\\n";
	    print SH "\t--start 1 --stop 500 --nBins 499 -dt NONE \\\n";
	    print SH "\t-R $bwaIndex{$reference{$project}} \\\n";
	    print SH "\t-I $outputDirectory/$project/bam/$sample.recal.bam \\\n";
	    print SH "\t-o $outputDirectory/$project/bam/$sample.recal.coverage\n\n";

	    print SH "# Identify areas for local realignment\n";
	    print SH "java -Xmx2G -jar /projects/p20742/tools/bin/GATK_v3.6//GenomeAnalysisTK.jar -T RealignerTargetCreator \\\n";
	    print SH "\t-R $bwaIndex{$reference{$project}} \\\n";
	    print SH "\t-I $outputDirectory/$project/bam/$sample.recal.bam \\\n";
	    print SH "\t-o $outputDirectory/$project/bam/$sample.realign.intervals \n\n";

	    print SH "# Realign problematic areas.\n";
	    print SH "java -Xmx2G -jar /projects/p20742/tools/bin/GATK_v3.6//GenomeAnalysisTK.jar -T IndelRealigner \\\n";
	    print SH "\t-R $bwaIndex{$reference{$project}} \\\n";
	    print SH "\t-known  $knownIndelsites{$reference{$sample}} \\\n";
	    print SH "\t-I $outputDirectory/$project/bam/$sample.recal.bam \\\n";
	    print SH "\t-o $outputDirectory/$project/bam/$sample.realigned.bam \\\n";
	    print SH "\t-targetIntervals $outputDirectory/$project/bam/$sample.realign.intervals \n\n";
	    
	 
	    print SH "# Replace uncorrected bam with corrected bam.\n";
	    print SH "cp $outputDirectory/$project/bam/$sample.bam  $outputDirectory/$project/bam/$sample.raw.bam\n";
	    print SH "cp $outputDirectory/$project/bam/$sample.bam.bai  $outputDirectory/$project/bam/$sample.raw.bam.bai\n"; 
 	    print SH "cp $outputDirectory\/$project\/bam\/$sample.realigned.bam $outputDirectory\/$project\/bam\/$sample.bam\n";

	    print SH "\n# Re-index sample.bam file to make sure the index is up-to-date.\n";
	    print SH "samtools index $outputDirectory\/$project\/bam\/$sample.bam\n";
 	    print SH "date\n\n";
	    
	    print SH "\n# Check if bamlist file exists, and if not, create it.\n";
	    print SH "if ! [ -f \"$outputDirectory\/$project\/bam\/bamlist.txt\" ];\nthen\necho \"$outputDirectory\/$project\/bam\/$sample.bam\" | cat > $outputDirectory\/$project\/bam\/bamlist.txt \n";
	    print SH "else\necho \"$outputDirectory\/$project\/bam\/$sample.bam\" | cat >> $outputDirectory\/$project\/bam\/bamlist.txt \nfi\n";
 	    if ($type eq "chipseq"){
 		if ($makeTracks == 1){
 		    print SH "# Make ChIPseq tracks.\n";
 		    # The multi mapping argument is not used right now.
 		    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
		    my $assembly = $reference{$sample};
		    if ($assembly eq "hg38.mp"){$assembly = "hg38";}
 		    print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$assembly --bamDir=$bamDirectory --sample=$sample --extLen=150\n";
 		    print SH "date\n\n";
 		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
 		    if ($uploadASHtracks == 1){
 			print SH "# Move tracks into Shilatifard directory structure\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
 			print SH "cp $outputDirectory\/$project\/bam\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			if ($uploadBAM == 1){
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "# Note that these files are not visible to browser unless you \"ake public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/ash-tracks/TANGO/$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://$s3path/$scientist.$project/ --region us-west-2\n";
			}
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
 		    }
		    if ($uploadPulmtracks == 1){
			if ($uploadBAM == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/m-328-data$scientist/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://m-328-data/$scientist.$project/ --region us-west-2\n";
			}
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp $outputDirectory\/$project\/bam\/$sample.bw s3://m-328-data/$scientist.$project/ --region us-west-2\n";
 		    }
		    if ($uploadLutracks == 1){
			if ($uploadBAM == 1){
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			    print SH "\n# Copy bamfiles to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    print SH "# To load, paste the following url into the custom tracks field:\n";
			    print SH "# http://s3-us-west-2.amazonaws.com/lwwlab/$scientist.$project/$sample.bam\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam s3://lwwlab/$scientist.$project/ --region us-west-2\n";
			    print SH "aws s3 cp $outputDirectory/$project/bam/$sample.bam.bai s3://lwwlab/$scientist.$project/ --region us-west-2\n";
			}
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
 			print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
 			print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
 			print SH "aws s3 cp $outputDirectory\/$project\/bam\/$sample.bw s3://lwwlab/$scientist.$project/ --region us-west-2\n";
 		    }
 		    print SH "mv $outputDirectory\/$project\/bam\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
 		    print SH "\n# Check if bwlist file exists, and if not, create it.\n";
 		    print SH "if  ! [-f \"$outputDirectory\/$project\/tracks\/bwlist.txt\"];\nthen\necho \"$outputDirectory\/$project\/tracks\/$sample.bw\" | cat > $outputDirectory\/$project\/tracks\/bwlist.txt \n";
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
			    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
#			    $moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
#			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
 			    print SH "samtools view -h $bamDirectory\/$sample.bam | awk '!(\$3 == \"$chr\" && \$4 > $start && \$4 < $stop){print \$0}' | samtools view -Sb - > $bamDirectory\/$sample.noVP.bam\n";
 			    print SH "\n# Make 4C tracks.\n";
 			    # The multi mapping argument is not used right now.
 			    # This is because Bowtie doesn't fill in the NH tag in the BAM file.
			    my $assembly = $reference{$sample};
			    if ($assembly eq "hg38.mp"){$assembly = "hg38";}
 			    print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$assembly --bamDir=$bamDirectory --sample=$sample.noVP --extLen=0\n";
 			    print SH "mv $bamDirectory\/$sample.noVP.bw $bamDirectory\/$sample.bw\n";
 			    print SH "date\n";
 			}
 		    } else {
 			print SH "\n# Make 4C tracks.\n";
 			# The multi mapping argument is not used right now.
 			# This is because Bowtie doesn't fill in the NH tag in the BAM file.
 			print SH "Rscript $NGSbartom/tools/bin/createChIPtracks.R --assembly=$reference{$sample} --bamDir=$bamDirectory --sample=$sample --extLen=0\n";
 		    }
 		    print SH "date\n\n";
 		    print SH "mkdir $outputDirectory\/$project\/tracks\n";
 		    if ($uploadASHtracks == 1){
 			print SH "\n# Move tracks into Shilatifard directory structure\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist\n";
 			print SH "mkdir /projects/b1025/tracks/TANGO/$scientist/$scientist.$project\n";
 			print SH "cp $outputDirectory\/$project\/bam\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project/\n";
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
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
	 if ($scheduler eq "MOAB"){
	     print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec msub {} ./ \\\;\n";
	 } elsif ($scheduler eq "SLURM"){
	     print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec msub {} ./ \\\;\n";
	 }
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
	my $result = "";
	if ($scheduler eq "MOAB"){
	    $result = `find $outputDirectory/$project/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec msub {} ./ \\\;`;
	    $result =~ s/\s+/\:/g;
	    $result =~ s/^://;
	    $result =~ s/:$//;
	} elsif ($scheduler eq "SLURM"){
	    $result = `find $outputDirectory/$project/scripts/ -iname \"run\_*\_$aligner\_align.sh\" -exec sbatch {} ./ \\\;`;
	    $result =~ s/Submitted batch job //g;
	    $result =~ s/\s+/\:/g;
	    $result =~ s/^://;
	    $result =~ s/:$//;
	    print STDERR "Predicted SLURM job IDs:  $result\n";
	}
	&datePrint( "Need to wait for the following jobs to finish:\n$result");
	open (SH, ">$outputDirectory/$project/scripts/AlignmentDependentScript.sh");
	print SH $header;
	my (%modulesLoaded, $moduleText);
	if ($scheduler eq "MOAB"){
	    print SH "#MSUB -N PostAlignmentAnalysis\n";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    print SH "#MSUB -W depend=afterok:$result\n";
	} elsif ($scheduler eq "SLURM"){
	    print SH "#SBATCH --job-name=PostAlignmentAnalysis\n";
	    print SH "#SBATCH --nodes=1\n";
	    print SH "#SBATCH -n $numProcessors\n";
	    print SH "#SBATCH --dependency=afterok:$result\n";
	}
	print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	print SH "module load deeptools/3.1.1\n";
	print SH "\necho \"Alignment jobs $result have finished.\"\n";
	my @alignJobs = split(/\:/,$result);
	if ($aligner eq "tophat"){
	    print SH "module load R/3.2.2\n";
	    print VER "EXEC module load R/3.2.2\n";
	    $moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
	    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	 #   $moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	 #   if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "\n# Make Tophat report summarizing alignment.\n";
	    print SH "Rscript $NGSbartom/tools/bin/createTophatReport.R --topHatDir=$outputDirectory\/$project\/Tophat_aln --nClus=$numProcessors\n";
	}elsif ($aligner eq "bowtie"){
	    print SH "find $bamDirectory\/*bowtie.log -type f -print -exec cat {} \\\; >> $outputDirectory\/$project\/alignlog.txt\n";
	}
	if ($aligner eq "star"){
	    print SH "touch $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "for r in $outputDirectory\/$project\/STAR_aln\/*.final.out\n";
	    print SH "do\n";
	    print SH "\techo \$r | cat >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Number of input reads\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Uniquely mapped reads number\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Uniquely mapped reads \%\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Number of reads mapped to multiple loci\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"\% of reads mapped to multiple loci\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Number of reads mapped to too many loci\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"\% of reads mapped to too many loci\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"\% of reads unmapped:\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"Number of chimeric reads\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\tgrep \"% of chimeric reads\" \$r >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "\techo \"=========================================================================\" | cat >> $outputDirectory\/$project\/alignlog.txt\n";
	    print SH "done\n";
	}
	if ($type eq "chipseq"){
	    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
	    print SH "\n# Plot ChIP fingerprint\n";
	    print SH "if ! [ -f $outputDirectory\/$project\/scripts\/plot_fingerprint.sh ];\nthen\npython3 /projects\/p20742\/tools\/bin\/getFingerprint.py -i $outputDirectory\/$project\/bam\/ -o $outputDirectory\/$project\/scripts\/\nsh $outputDirectory\/$project\/scripts\/plot_fingerprint.sh\n";
	    #    print SH "wait\nmv $outputDirectory\/$project\/scripts\/fingerprint.pdf $outputDirectory\/$project\/bam\/\nfi\n";
	}
	if ($scheduler eq "SLURM"){
	    print SH "\n # Use Seff to report on job outcome.\n";
	    foreach my $job (@alignJobs){
		print SH "seff $job > $outputDirectory/metadata/SLURM.align.$job.seff.txt\n";
	    }
	}
	close SH;
	&datePrint("Creating dependent job that will only run after alignments finish.");
	my $result2 = "";
	if ($scheduler eq "MOAB"){
	    $result2 = `msub $outputDirectory/$project/scripts/AlignmentDependentScript.sh`;
	    $result2 =~ s/\s+//g;
	} elsif ($scheduler eq "SLURM"){
	    $result2 = `sbatch $outputDirectory/$project/scripts/AlignmentDependentScript.sh`;
	    $result2 =~ s/Submitted batch job //g;
	    $result2 =~ s/\s//g;
#	    print STDERR "Predicted SLURM job ID:  $result2\n";
	}
	my $jobfinished = "no";
	# Wait until the job is Complete.
	&datePrint("Waiting for job $result2 to finish. (each . = 300 seconds)");
	waitForJob($result2);
	print STDERR "\n";
	&datePrint("Job $result2 done.  Continuing.");
    }
}

if ($buildSampleCheck == 1){
    &datePrint("Creating scripts for checking sample identity.");
    my @samples;
    # Foreach project (TANGO/MOLNG):
    foreach my $project (keys(%samples)){
	my $cmd = "mkdir $outputDirectory\/$project\/SampleID";
	system($cmd);
	my $shScript = "$outputDirectory\/$project\/scripts\/run\_$project\_SampleCheck.sh";
	&datePrint("Printing to $shScript");
	open (SH,">$shScript");
	my (%modulesLoaded, $moduleText);
	print SH "$header";
	if ($scheduler eq "MOAB"){
	    print SH "#MSUB -N $project\_sampleCheck\n";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	} elsif ($scheduler eq "SLURM"){
	    print SH "#SBATCH --job-name=$project\_sampleCheck\n";
	    print SH "#SBATCH --nodes=1\n";
	    print SH "#SBATCH -n $numProcessors\n";
	}
	print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n";
	print SH "\n# Setup paths for the analysis\n";
	print SH "export PATH=\$PATH:$NGSbartom/tools/bin/\n";
	print SH "export NCM_HOME=$NGSbartom/tools/bin/NGSCheckMate/\n\n";
	$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
	if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
	my $fastqList = "$outputDirectory\/$project\/SampleID\/$project.fastqList.txt";
	open (FQL,">$fastqList");
	@samples = uniq(split(/\,/,$samples{$project}));
	# Foreach sample within the project:
	foreach my $sample (@samples){
	    my @fastqs = split(/\,/,$fastqs{$sample});
	    my $cmd = "";
	    if (!(-e "$outputDirectory\/$project\/fastq")){
		$cmd .= "mkdir $outputDirectory\/$project\/fastq\n";
	    }
	    my $oldfastqs = "@fastqs";
	    my @newfastqs = "";
	    for my $fastq (@fastqs){
		my $basefilename = basename($fastq);
		$fastq = "$outputDirectory\/$project\/fastq\/$basefilename";
		push (@newfastqs, $fastq);
	    }
	    my $newfastqs = "@newfastqs";
	    #	    $newfastqs =~ s/ /,/g;
	    if ($oldfastqs ne $newfastqs){
		print STDERR "Old fastqs: $oldfastqs\nNew fastqs: $newfastqs\n";
		system($cmd);
	    }
	    @fastqs = split(/\s+/,$newfastqs);
	    for my $fastq (@fastqs){
		my $fastqLabel = basename($fastq);
		$fastqLabel =~ s/\.fastq.gz$//g;
		$fastqLabel =~ s/\.fq.gz$//g;
		$fastqLabel =~ s/\.fastq$//g;
		$fastqLabel =~ s/\.fq$//g;
		if ($fastq ne ""){
		    print FQL "$fastq\t$fastqLabel\n";
		}
	    }
	    print SH "\n# Run NGSCheckmate on the fastq files listed in  $outputDirectory\/$project\/SampleID\/$project.fastqList.txt\n";
	    print SH "python \$NCM_HOME/ncm_fastq.py -l $outputDirectory\/$project\/SampleID\/$project.fastqList.txt -pt \$NCM_HOME/SNP/SNP.pt -O $outputDirectory\/$project\/SampleID/NGSCheckmateResults >& $outputDirectory\/$project\/SampleID\/$project.ncm.fq.log\n\n";
	    close(SH);
	}
	if ($runSampleCheck == 0){
	    # Print tips on running the tophat shell scripts.
	    print STDERR "To execute all scripts, use the following command:\n";
	    # "$outputDirectory\/$project\/scripts\/run\_$project\_SampleCheck.sh";
	    if ($scheduler eq "MOAB"){
		print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_$project\_SampleCheck.sh\" -exec msub {} ./ \\\;\n";
	    } elsif ($scheduler eq "SLURM"){
		print STDERR "find $outputDirectory/*/scripts/ -iname \"run\_$project\_SampleCheck.sh\" -exec sbatch {} ./ \\\;\n";
	    }
	} elsif ($runSampleCheck == 1){
	    &datePrint("Submitting job for SampleCheck.");
	    if ($scheduler eq "MOAB"){
		$cmd = "msub $outputDirectory/$project/scripts/run\_$project\_SampleCheck.sh";
	    } elsif ($scheduler eq "SLURM"){
		$cmd = "sbatch $outputDirectory/$project/scripts/run\_$project\_SampleCheck.sh";
	    }
	    system($cmd);
	}
    }
}

if ($runRNAstats == 1){
    if ($type eq "RNA"){
	foreach my $project (keys(%samples)){
	    open (SH, ">$outputDirectory/$project/scripts/runRNAstats_summary.sh");
	    my (%modulesLoaded, $moduleText);
	    print SH $header;
	    if ($scheduler eq "MOAB"){
		print SH "#MSUB -N $project\_runRNAstats\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		print SH "#MSUB -l walltime=48:00:00\n";
	    } elsif ($scheduler eq "SLURM"){
		print SH "#SBATCH --job-name=$project\_runRNAstats\n";
		print SH "#SBATCH --nodes=1\n";
		print SH "#SBATCH -n $numProcessors\n";
		print SH "#SBATCH --time=48:00:00\n";
	    }
	    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	    print SH "module load R/3.2.2\n";
	    print VER "EXEC module load R/3.2.2\n";
	    $moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
	    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
	    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	    #$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	    #if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
	    print SH "module unload mpi/openmpi-1.6.3-gcc-4.6.3\n";
	    print VER "EXEC module unload mpi/openmpi-1.6.3-gcc-4.6.3\n";
	    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
	    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
	    print SH "module load parallel\n";
	    print VER "EXEC RSeQC\n";
	    my $strandrule = "";
	    if ($stranded == 0){
		$strandrule = "none";
	    } else {
		if ($runPairedEnd == 1){
		    $strandrule = "1++,1–,2+-,2-+";
		} else {
		    $strandrule = "+-,-+";
		}
	    }
	    my @samples = uniq(split(/\,/,$samples{$project}));
	    foreach my $sample (@samples){
		open (SSH, ">$outputDirectory/$project/scripts/runRNAstats.$sample.sh");
		print SSH $header;
		my (%modulesLoaded, $moduleText);
		if ($scheduler eq "MOAB"){
		    print SSH "#MSUB -N $project\_runRNAstats.$sample\n";
		    print SSH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		} elsif ($scheduler eq "SLURM"){
		    print SSH "#SBATCH --job-name=$project\_runRNAstats.$sample\n";
		    print SSH "#SBATCH --nodes=1\n";
		    print SSH "#SBATCH -n $numProcessors\n";
		}
		print SSH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
		print SSH "module load R/3.2.2\n";
		$moduleText = &checkLoad("bowtie2/2.2.6",\%modulesLoaded);
		if ($moduleText ne ""){	print SSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bowtie2/2.2.6"} = 1;}
		$moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		if ($moduleText ne ""){ print SSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		#$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
		#if ($moduleText ne ""){ print SSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print SSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
		print SSH "module unload mpi/openmpi-1.6.3-gcc-4.6.3\n";
		print SSH "module load parallel\n";
		if (!(-e "$outputDirectory\/$project\/bam\/$sample.bam")){
		    print SSH "\n# When this script was generated, $sample.bam did not exist in the bam directory.  If that is still the case when this script is run, it will fail.  In that case, please align and re-run.\n";
		}
		print SSH "\n# Use RSeQC to infer experiment type for $sample.\n";
		print SSH "infer_experiment.py -r $samplegenebed{$reference{$project}} -i $outputDirectory\/$project\/bam\/$sample.bam > $outputDirectory\/$project\/bam\/$sample\_inferredExperiment.txt\n";
		print SSH "\n# Use RSeQC to check saturation for $sample.\n";
		print SSH "RPKM_saturation.py -r $samplegenebed{$reference{$project}} -i $outputDirectory\/$project\/bam\/$sample.bam -o $outputDirectory\/$project\/bam\/$sample\n";
		print SSH "\n# Examine read duplication for $sample.\n";
		print SSH "read_duplication.py -i $outputDirectory\/$project\/bam\/$sample.bam -o $outputDirectory\/$project\/bam\/$sample\n";
		print SSH "\n# How many novel splice sites are found in the data? (Only relevant for STAR aligned data)\n";
		print SSH "\njunction_annotation.py -i $outputDirectory\/$project\/bam\/$sample.bam -r $samplegenebed{$reference{$project}} -o $outputDirectory\/$project\/bam\/$sample\n";
		print SSH "\n# Check splice site saturation for $sample.\n";
		print SSH "\njunction_saturation.py -i $outputDirectory\/$project\/bam\/$sample.bam -r $samplegenebed{$reference{$project}} -o $outputDirectory\/$project\/bam\/$sample\n";
		print SSH "\n# If sorted bam doesn't exist, create it.\n";
		print SSH "if ! [ -f $outputDirectory\/$project\/bam\/$sample.sorted.bam\" ]; then\n";
		print SSH "\tsamtools sort $outputDirectory\/$project\/bam\/$sample.bam > $outputDirectory\/$project\/bam\/$sample.sorted.bam\n";
		print SSH "\tsamtools index $outputDirectory\/$project\/bam\/$sample.sorted.bam\n";
		print SSH "fi\n";
		print SSH "\n# Estimate TIN (transcript integrity number) for each transcript.\n";
		print SSH "tin.py -i $outputDirectory\/$project\/bam\/$sample.sorted.bam -r $samplegenebed{$reference{$project}} \n";
		close(SSH);
		my $result2 = "";
		if ($scheduler eq "MOAB"){
		    $result2 = `msub $outputDirectory/$project/scripts/runRNAstats.$sample.sh`;
		} elsif ($scheduler eq "SLURM"){
		    $result2 = `sbatch $outputDirectory/$project/scripts/runRNAstats.$sample.sh`;
		    $result2 =~ s/Submitted batch job //g;
		}
	    }
	    print SH "\n# Use RSeQC to estimate gene body coverage across all samples.\n";
	    print SH "geneBody_coverage.py -r $samplegenebed{$reference{$project}} -i $outputDirectory\/$project\/bam\/ -o $outputDirectory\/$project\/bam\/$project\n";
	    &datePrint("Launching RNA stats script for project $project");
	    my $result3 = "";
	    if ($scheduler eq "MOAB"){
		$result3 = `msub $outputDirectory/$project/scripts/runRNAstats_summary.sh`;
	    } elsif ($scheduler eq "SLURM"){
		$result3 = `sbatch $outputDirectory/$project/scripts/runRNAstats_summary.sh`;
	    }
	}
    } else {
	&datePrint("Skipping runRNAstats for nonRNA project");
    }
}

my %exomePairing;
if ($buildGenotyping ==1) {
    if (($type eq "RNA") || ($type eq "exome")){
	foreach my $project (keys(%samples)){
	    my @samples = uniq(split(/\,/,$samples{$project}));
	    if (-d "$outputDirectory\/$project\/genotype/") {
		print STDERR "Genotype Directory found.\n";
	    }else {system("mkdir $outputDirectory\/$project\/genotype");print STDERR "Genotype directory created.\n";}
	    if ($exomeDescription ne ""){
		open(IN,$exomeDescription);
		print STDERR "Found Exome pairing file $exomeDescription\n";
	        while(<IN>){
		    chomp $_;
		    my($exome,$normal) = split(/\s+/,$_);
		    $exome =~ s/.bam$//g;
		    $normal =~ s/.bam$//g;
#		    print STDERR "Pairing $exome and $normal\n";
		    $exomePairing{$exome} = $normal;
		}
	    }
	    # Foreach sample within the project:
	    foreach my $sample (@samples){
		open (SH, ">$outputDirectory/$project/scripts/run_$sample\_genotype.sh");
		my (%modulesLoaded, $moduleText);
		print SH $header;
		if ($scheduler eq "MOAB"){
		    print SH "#MSUB -N $sample\_genotype\n";
		    #print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		    # Some of the steps below can be parallelized, but not all,
		    # and more experimentation is required for optimization.  In
		    # the meantime, I'm setting a flat number of processors of 6
		    # for genotyping purposes.
		    # 2019-09-07; removing the flat number of processors.
#		    print SH "#MSUB -l nodes=1:ppn=6\n";
		    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		} elsif ($scheduler eq "SLURM"){
		    print SH "#SBATCH --job-name=$sample\_genotype\n";
		    print SH "#SBATCH --nodes=1\n";
#		    print SH "#SBATCH -n 6\n";
		    print SH "#SBATCH -n $numProcessors\n";
		}
		if ($type eq "RNA"){
		    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
		    my $PICARD = "/software/picard/1.131/picard-tools-1.131/picard.jar";
		    print SH "\n\n";
		    print SH "# Sort BAM file.\n";
		    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		    print SH "samtools sort -o $outputDirectory\/$project\/bam\/$sample.sorted.bam $outputDirectory\/$project\/bam\/$sample.bam\n";
		    print SH "date\n\n";
		    $moduleText = &checkLoad("picard/1.131",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"picard/1.131"} = 1;}
		    print SH "# Mark Duplicates with Picard.\n";
		    print SH "java -jar $PICARD MarkDuplicates \\\n";
		    print SH "\tI=$outputDirectory\/$project\/bam\/$sample.sorted.bam \\\n";
		    print SH "\tO=$outputDirectory\/$project\/bam\/$sample.mdup.bam \\\n";
		    print SH "\tM=$outputDirectory\/$project\/bam\/$sample.mdup.metrics.txt\n";
		    print SH "date\n\n";
		    print SH "# Reorder mdup BAM file with Picard.\n";
		    print SH "java -jar $PICARD ReorderSam \\\n";
		    print SH "\tI=$outputDirectory\/$project\/bam\/$sample.mdup.bam \\\n";
		    print SH "\tO=$outputDirectory\/$project\/bam\/$sample.mdup.reordered.bam \\\n";
		    print SH "\tR=$gatkRef{$reference{$sample}} \\\n";
		    print SH "\tCREATE_INDEX=true\n";
		    print SH "date\n\n";
		    print SH "# GATK 4 is not working.  Reverting to GATK 3.7\n";
		    #		print SH "# Set up path for GATK 3.7\n";
		    $moduleText = &checkLoad("gatk/3.7.0",\%modulesLoaded);
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"gatk/3.7.0"} = 1;}
		    print SH "export PATH=\$PATH:/software/gatk/3.7.0/\n\n";
		    print SH "# Split Reads at splicing events (runs of Ns in CIGAR string)\n";
		    #		print SH "# Note that in earlier versions of gatk / aligners there was a problem with low quality scores \n";
		    #		print SH "# at splice junctions but that is now handled automatically by gatk4.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T SplitNCigarReads \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I $outputDirectory\/$project\/bam\/$sample.mdup.reordered.bam \\\n";
		    print SH "\t-o $outputDirectory\/$project\/bam\/$sample.split.bam \\\n";
		    print SH "\t-U ALLOW_N_CIGAR_READS -fixNDN\n";
		    print SH "date\n\n";
		    print SH "# Find Target regions for Realignment.\n";  # Not available in gatk4, only gatk3.
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T RealignerTargetCreator \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I $outputDirectory\/$project\/bam\/$sample.split.bam \\\n";
		    print SH "\t-o $outputDirectory\/$project\/bam\/$sample.split.intervals.list \\\n";
		    print SH "\t--known $knownIndelsites{$reference{$sample}}\n";
		    print SH "date\n\n";
		    print SH "# Realign indels in target regions.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T IndelRealigner \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I $outputDirectory\/$project\/bam\/$sample.split.bam \\\n";
		    print SH "\t-targetIntervals $outputDirectory\/$project\/bam\/$sample.split.intervals.list \\\n";
		    print SH "\t-known $knownIndelsites{$reference{$sample}} \\\n";
		    print SH "\t-o $outputDirectory\/$project\/bam\/$sample.split.real.bam\n";
		    print SH "date\n\n";
		    print SH "# Generating Base Recalibration Table.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T BaseRecalibrator \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I $outputDirectory\/$project\/bam\/$sample.split.real.bam \\\n";
		    print SH "\t-o $outputDirectory\/$project\/genotype\/$sample.split.real.recal.table \\\n";
		    print SH "\t-knownSites $knownSNPsites{$reference{$sample}}\n";
		    print SH "date\n\n";
		    #		print SH "##Commenting out HaplotypeCaller as Mutect2 is working better.\n";
		    print SH "# Calling SNPs and Indels with HaplotypeCaller.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T HaplotypeCaller \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I $outputDirectory\/$project\/bam\/$sample.split.real.bam \\\n";
		    print SH "\t-o $outputDirectory\/$project\/genotype\/$sample.raw.snps.indels.vcf \\\n";
		    print SH "\t--dbsnp $knownSNPsites{$reference{$sample}}\n";
		    print SH "date\n\n";
		    print SH "# Filter HC calls.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T SelectVariants \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t--excludeFiltered \\\n";
		    print SH "\t--variant $outputDirectory\/$project\/genotype\/$sample.raw.snps.indels.vcf \\\n";
		    print SH "\t--out $outputDirectory\/$project\/genotype\/$sample.filtered.snps.indels.hc.vcf\n";
		    print SH "date\n\n";
		    print SH "# Calling SNPs and Indels with Mutect2.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T MuTect2 \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t-I:tumor $outputDirectory\/$project\/bam\/$sample.split.bam \\\n";
		    print SH "\t--out $outputDirectory\/$project\/genotype\/$sample.raw.snps.indels.m2.vcf \\\n";
		    print SH "\t--dbsnp $knownSNPsites{$reference{$sample}}\n";
		    print SH "date\n\n";
		    print SH "# Filter Mutect2 calls.\n";
		    print SH "java -jar /software/gatk/3.7.0/GenomeAnalysisTK.jar \\\n";
		    print SH "\t-T SelectVariants \\\n";
		    print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
		    print SH "\t--excludeFiltered \\\n";
		    print SH "\t--variant $outputDirectory\/$project\/genotype\/$sample.raw.snps.indels.m2.vcf \\\n";
		    print SH "\t--out $outputDirectory\/$project\/genotype\/$sample.filtered.snps.indels.m2.vcf\n";
		    print SH "date\n\n";
#		    print SH "# Sort and remove duplicates from VCF file.\n";
#		    print SH "sort $outputDirectory\/$project\/genotype\/$sample.filtered.snps.indels.m2.vcf | uniq > $outputDirectory\/$project\/genotype\/$sample.filtered.snps.indels.m2.sorted.vcf\n";
#		    close SH;
		} elsif ($type eq "exome"){
		    ## Return Here.
		    print SH "# Index BAM file.\n";
		    $moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		    if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"gatk/3.6.0"} = 1;}
		    print SH "export PATH=\$PATH:/software/gatk/3.6.0/\n\n";
#		    print SH "samtools index -b $outputDirectory\/$project\/bam\/$sample.bam\n";
#		    print SH "date\n\n";
		    my $normal = "";
		    if (exists($exomePairing{$sample})){
			$normal = $exomePairing{$sample};
#			print SH "samtools index -b $outputDirectory\/$project\/bam\/$normal.bam\n";
#			print SH "date\n\n";
			print SH "# Call variants with MuTect2\n";
			print SH "java -Xmx2G -jar /software/gatk/3.6.0/GenomeAnalysisTK.jar -T MuTect2 \\\n";
			print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
			print SH "\t-I:tumor $outputDirectory\/$project\/bam\/$sample.bam \\\n";
			print SH "\t-I:normal $outputDirectory\/$project\/bam\/$normal.bam \\\n";
			print SH "\t--dbsnp $knownSNPsites{$reference{$sample}}\\\n";
			print SH "\t--cosmic $cosmicSNPsites{$reference{$sample}}\\\n";
			print SH "\t--output_mode EMIT_VARIANTS_ONLY \\\n";
			print SH "\t-o $outputDirectory\/$project\/genotype\/$sample.relativeTo$normal.vcf.gz\\\n\n";
			print SH "# Filter variants\n";
			print SH "java -Xmx2G -jar /software/gatk/3.6.0/GenomeAnalysisTK.jar -T VariantFiltration \\\n";
			print SH "\t-V $outputDirectory\/$project\/genotype\/$sample.relativeTo$normal.vcf.gz\\\n";
			print SH "\t-O $outputDirectory\/$project\/genotype\/$sample.relativeTo$normal.filtered.vcf.gz\\\n";		
		    } else { 
			print SH "# No normal control found for sample $sample.\n";
			print SH "# Call variants with MuTect2\n";
			print SH "java -Xmx2G -jar /software/gatk/3.6.0/GenomeAnalysisTK.jar -T MuTect2 \\\n";
			print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
			print SH "\t-I:tumor $outputDirectory\/$project\/bam\/$sample.bam \\\n";
			print SH "\t--dbsnp $knownSNPsites{$reference{$sample}}\\\n";
			print SH "\t--cosmic $cosmicSNPsites{$reference{$sample}}\\\n";
			print SH "\t--output_mode EMIT_VARIANTS_ONLY \\\n";
			print SH "\t-o $outputDirectory\/$project\/genotype\/$sample.vcf.gz\\\n\n";
			print SH "# Filter variants\n";
			print SH "java -Xmx2G -jar /software/gatk/3.6.0/GenomeAnalysisTK.jar -T VariantFiltration \\\n";
			print SH "\t-V $outputDirectory\/$project\/genotype\/$sample.vcf.gz\\\n";
			print SH "\t-R $gatkRef{$reference{$sample}} \\\n";
			print SH "\t-o $outputDirectory\/$project\/genotype\/$sample.filtered.vcf.gz\\\n";		
		    }
		} else { die "ERR: Genotyping does not appear to be implemented for type $type.\n";}
	    }
	    my $result = "";
	    if ($scheduler eq "MOAB"){
		if ($runGenotyping == 1){
		    &datePrint("Submitting jobs for genotype analysis.");
		    $result = `find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec msub {} ./ \\\;`;
		} else {
		    print STDERR "Not submitting genotype scripts.  To run them use: \`find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec msub \{\} ./ \\\;\`\n";
		}
	    } elsif ($scheduler eq "SLURM"){
		if ($runGenotyping == 1){
		    &datePrint("Submitting jobs for genotype analysis.");
		    $result = `find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec sbatch {} ./ \\\;`;
		} else {
		    print STDERR "Not submitting genotype scripts.  To run them use: \`find $outputDirectory/$project/scripts/ -iname \"*genotype.sh\" -exec sbatch \{\} ./ \\\;\`\n";
		}
	    }
#    } else { print "ERR: Genotyping not yet implemented for whole genome or exome sequencing.\n";}
	}
    }
}


if ($buildEdgeR ==1) {
    foreach my $project (keys(%samples)){
	my @countsfiles;
	open (SH, ">$outputDirectory/$project/scripts/downstreamRNAanalysis.sh");
	my (%modulesLoaded, $moduleText);
	print SH $header;
	if ($scheduler eq "MOAB"){
	    print SH "#MSUB -N downstreamRNAanalysis\n";
	    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	} elsif ($scheduler eq "SLURM"){
	    print SH "#SBATCH --job-name=downstreamRNAanalysis\n";
	    print SH "#SBATCH --nodes=1\n";
	    print SH "#SBATCH -n $numProcessors\n";
	}
	print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	print SH "module load R/3.2.2\n";
	if ($type eq "RNA"){
	    if ($startFromBAM == 1){
		# If the bam directory was specified in the input (no base space directory or fastq directory)
		if ($htseq == 1 || $bedtools == 1) {
		    my @samples = uniq(split(/\,/,$samples{$project}));
		    # Foreach sample within the project:
		    if ($htseq == 1) {
		    	print SH "module unload mpi\n";
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
			print VER "EXEC htseq 0.6.1\n";
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
			if ($modulesLoaded{"python/anaconda"} == 1){ print SH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
	    		print SH "module load gcc/4.8.3\n";
		    }
		    if ($bedtools == 1) {
			$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
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
			$modulesLoaded{"bedtools/2.17.0"} = 0;
		    }
		}
		print SH "module load R/3.2.2\n";
		print VER "EXEC module load R/3.2.2\n";
		if ($makeTracks == 1){
		    my @samples = uniq(split(/\,/,$samples{$project}));
		    # Foreach sample within the project:
		    foreach my $sample (@samples){
		    # if the user wants tracks, put that in the downstream analysis script.
			print SH "# Create RNA seq Tracks\n";
			my $assembly = $reference{$sample};
			if ($assembly eq "hg38.mp"){$assembly="hg38";}
			print SH "Rscript $NGSbartom/tools/bin/createRNAseqTracks3.R --assembly=$assembly --bamDir=$bamDirectory\/ --sample=$sample --stranded=$stranded  --multiMap=$multiMap\n";
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
			    if ($stranded == 1){
				print SH "cp $bamDirectory\/$sample.minus.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
				print SH "cp $bamDirectory\/$sample.plus.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			    } else {
				print SH "cp $bamDirectory\/$sample.bw /projects/b1025/tracks/TANGO/$scientist\/$scientist.$project\/\n";
			    }
			    $moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			    if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}

			    print SH "\n# Copy bigwigs to Amazon S3, for UCSC genome browser to access.\n";
			    print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
			    if ($stranded == 1){
				print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.minus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
				print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.plus.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    } else {
				print SH "aws s3 cp /projects/b1025/tracks/TANGO/$scientist/$scientist.$project/$sample.bw s3://$s3path/$scientist.$project/ --region us-west-2\n";
			    }
			} else {
			    if ($stranded == 1){
				print SH "mv $bamDirectory\/$sample.minus.bw $outputDirectory\/$project\/tracks\/\n";
				print SH "mv $bamDirectory\/$sample.plus.bw $outputDirectory\/$project\/tracks\/\n";
			    } else {
				print SH "mv $bamDirectory\/$sample.bw $outputDirectory\/$project\/tracks\/\n";
			    }
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
		my $assembly = $reference{$project};
		if ($assembly eq "hg38.mp"){$assembly="hg38";}
		print SH "Rscript $NGSbartom/tools/bin/createRNAcounts.R --assembly=$assembly --bamDir=$bamDirectory --numCores=$numProcessors --txdbfile=$txdbfile{$reference{$project}}\n";
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
	    	print SH "perl $NGSbartom/tools/bin/makeHTseqCountsTable.pl $bamDirectory\/ $gff{$reference{$project}} $bamDirectory\/\n";
		push (@countsfiles,"$bamDirectory/htseq.all.counts.txt");
	    }
	    if ($rsem == 1) {
		push (@methods,"rsem");
		print SH "\n# Make RSEM counts table.\n";
	    	print SH "perl $NGSbartom/tools/bin/makeRSEMcountsTable.pl $bamDirectory\/ $gff{$reference{$project}} $bamDirectory\/ isoforms\n";
		push (@countsfiles,"$bamDirectory/rsem.all.counts.txt");
	    }
	    if ($bedtools == 1) {
		push (@methods,"bedtools");
		print SH "\n# Make Bedtools counts table.\n";
	    	print SH "perl $NGSbartom/tools/bin/makeBEDtoolsCountsTable.pl $bamDirectory\/ $bamDirectory\/\n";
		push (@countsfiles,"$bamDirectory/htseq.all.counts.txt");
	    }
	    if ($comparisons ne ""){
		foreach my $method (@methods){
		    print SH "\n# Run EdgeR, using comparisons file, without MDS plot (which sometimes crashes), $method.\n";
		    if ($reference{$project} =~ /\.mp/){
			print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --comparisonFile=$comparisons --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis --runMDS=0 --filterOff=1\n";
		    } else {
			print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --comparisonFile=$comparisons --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis --runMDS=0\n";
		    }
		    print SH "\n# Run EdgeR, creating MDS plot, but not running comparisons, $method.\n";
		    print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis --runMDS=1\n";
		    print SH "date\n";
		    if ($method eq "granges"){
			print SH "\n# Create Correlation plot for all samples\n";
			print SH "# This plot is still in development, so don't over-interpret, $method.\n";
			print SH "Rscript $NGSbartom/tools/bin/makeCorrelationPlotAllSamples.R --countFile=$bamDirectory\/$method.all.counts.rda --outputDirectory=$outputDirectory\/$project\/analysis\n";
		    }
		}
		print SH "date\n";
		print SH "# Processing $comparisons file\n";
		open(CMP,$comparisons);
		while(<CMP>){
		    chomp $_;
		    my ($comp,@groups) = split(/\,/,$_);
		    if ("@groups" =~ /^[1\-0\s]+/){
			foreach my $method (@methods){
			    print SH "\n# Create Gene Lists for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/bin/fromEdgeRtoGeneList.R --edgeRfile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01\n";
			    print SH "\n# Create MA plot for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/bin/makeMAplot.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --labelTop=1\n";
			    print SH "date\n";
			    print SH "\n# Create Big Heatmaps for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/bin/makeBigHeatmap.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --countFile=$outputDirectory\/$project\/analysis\/$method.normCounts.txt\n";
			    print SH "date\n";
			    print SH "\n# Run GO analysis for comparison $comp, method $method\n";
			    print SH "Rscript $NGSbartom/tools/bin/runGOforDEG.R --degFile=$outputDirectory\/$project\/analysis\/$comp.$method.edgeR.txt --adjp=0.01 --assembly=$reference{$project}\n";
			    print SH "date\n";
			}
		    }
		}
		close(CMP);
	    } else {
		foreach my $method (@methods){
		    print SH "\n# Create MDS plot for samples, with count method $method.\n";
		    print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --countFile=$bamDirectory\/$method.all.counts.txt --numCores=$numProcessors --outputDirectory=$outputDirectory\/$project\/analysis\n";
		}
	    }
	    if ($runEdgeR == 1){
		&datePrint("Submitting job for downstream RNA-seq analysis.");
		if ($scheduler eq "MOAB"){
		    $cmd = "msub $outputDirectory/$project/scripts/downstreamRNAanalysis.sh";
		} elsif ($scheduler eq "SLURM"){
		    $cmd = "sbatch $outputDirectory/$project/scripts/downstreamRNAanalysis.sh";
		}
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
	    if ($buildNGSplot == 1){
		$cmd .= "\nmkdir $outputDirectory\/$project\/analysis";
		$cmd .= "\nmkdir $outputDirectory\/$project\/NGSplotPipeline";
	    }
	    system($cmd);
	    my $broadPeakCheck = `grep -i \"broad\" $chipDescription`;
	    if ($broadPeakCheck ne ""){
		&datePrint("Found some broad peaks to call. Creating SICER script.");
		open (BSH, ">$outputDirectory/$project/scripts/runSICER\_callPeaks.sh");
		my (%modulesLoaded, $moduleText);
		print BSH $header;
		if ($scheduler eq "MOAB"){
		    print BSH "#MSUB -N callSicerPeaks\n";
		    print BSH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		} elsif ($scheduler eq "SLURM"){
		    print BSH "#SBATCH --job-name=callSicerPeaks\n";
		    print BSH "#SBATCH --nodes=1\n";
		    print BSH "#SBATCH -n $numProcessors\n";
		}
		print BSH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
		print BSH "\nmodule unload R\n";
		print BSH "module unload mpi\n";
		$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
		if ($moduleText ne ""){ print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
		$moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		if ($moduleText ne ""){	print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
	#	$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
	#	if ($moduleText ne ""){	print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
		$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
		if ($moduleText ne ""){ print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
		print BSH "export PATH=$NGSbartom/tools/bin/SICER_V1.1/SICER/:\$PATH\n";
		print VER "EXEC $NGSbartom/tools/bin/SICER_V1.1/SICER/\n";
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
			my (%modulesLoaded, $moduleText);
			if ($scheduler eq "MOAB"){
			    print SH "#MSUB -N $ip\_NarrowPeaks\n";
			    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
			} elsif ($scheduler eq "SLURM"){
			    print SH "#SBATCH --job-name=$ip\_NarrowPeaks\n";
			    print SH "#SBATCH --nodes=1\n";
			    print SH "#SBATCH -n $numProcessors\n";
			}
			print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
			#print SH "module load R/3.2.2\n";
			print SH "export PATH=\$PATH:$NGSbartom/tools/bin/MACS-1.4.2/bin\n";
			print VER "EXEC $NGSbartom/tools/bin/MACS-1.4.2\n";
			print SH "export PYTHONPATH=$NGSbartom/tools/bin/MACS-1.4.2/lib/python2.6/site-packages:\$PYTHONPATH\n";
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
			    print SH "$NGSbartom/tools/bin/bedToBigBed $outputDirectory\/$project\/peaks\/$ip.macsPeaks.capped.bed $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb\n";
			    print SH "date\n";
			    if ($uploadASHtracks == 1){
				$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
				if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
				print SH "\n# Copy bigBed to Amazon S3, for UCSC genome browser to access.\n";
				print SH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
				print SH "aws s3 cp $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb s3://$s3path/$scientist.$project/ --region us-west-2\n";
				print SH "date\n";
				if ($modulesLoaded{"python/anaconda"} == 1){ print SH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
			    }
			    print SH "\n# Create track description for bigBed file.\n";
			    print SH "echo \"track type=bigBed name=$ip.macsPeaks description=\\\"MACS peaks in $ip relative to $input\\\" graphtype=bar maxHeightPixels=128:60:11 visibility=dense color=0,0,0 itemRGB=on useScore=1 autoScale=on bigDataUrl=https://s3-us-west-2.amazonaws.com/$s3path/$scientist.$project/$ip.macsPeaks.bb\" | cat > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bb.header.txt\n";
			    print SH "date\n";
			}
			$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
			if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
			print SH "\n# Calculate coverage of reads at peaks.\n";
			print SH "bedtools coverage -b $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed -abam $outputDirectory\/$project\/bam\/$ip.bam >  $outputDirectory\/$project\/peaks\/$ip.macsPeaks.cov.txt\n";
			print SH "awk \'{printf \"\%s\\t\%d\\t%d\\t\%s\\t\%f\\n\", \$1,\$2,\$3\,\$4,\$5\/\$7\}\'  $outputDirectory\/$project\/peaks\/$ip.macsPeaks.cov.txt >  $outputDirectory\/$project\/peaks\/$ip.macsPeaks.cpm.bed\n";
			print SH "\n# Annotate peaks with nearby genes.\n";
			print VER "EXEC module load R/3.2.2\n";
			print SH "module load R/3.2.2\n";
			print SH "grep -v random $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.nonRandomChr.bed\n";
			print SH "Rscript $NGSbartom/tools/bin/addGenesToBed.R --peakFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.nonRandomChr.bed --outputDirectory=$outputDirectory\/$project\/peaks --assembly=$reference{$project} --txdbfile=$txdbfile{$reference{$project}}\n";
			print SH "\n# Center peaks and sort by peak width, from large to small.\n";
			print SH "perl $NGSbartom/tools/bin/NGSplotPipeline/NGSplotFilesScripts/convertToNGSplotSortedCenteredBED.pl $outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed $outputDirectory\/$project\/peaks\/$ip.sorted.centered.bed\n";
			print SH "\n# Build NGS plot config file.\n";
			if (-e "$outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultConfig.txt") {
			    print SH "rm $outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultConfig.txt\n";
			}
			print SH "touch $outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultConfig.txt\n";
			print SH "if ! [ -f \"$outputDirectory\/$project\/bam\/bamlist.txt\" ];\nthen\nls $outputDirectory\/$project\/bam\/*.bam > $outputDirectory\/$project\/bam\/bamlist.txt \nfi\n";
			print SH "date\n\n";
			print SH "for bam in \`cat $outputDirectory/$project/bam/bamlist.txt\`\n";
			print SH "do\n";
			print SH "\techo \$bam\n";
			print SH "\tbamsample=\$\(basename \$bam\)\n";
			print SH "\tprintf \"\%s\\t\%s\\t\%s\\n\" \"\$bam\" \"$outputDirectory\/$project\/peaks\/$ip.sorted.centered.bed\" \"\$bamsample\" \>\> $outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultConfig.txt\n";
			print SH "done\n";
			print SH "\nmodule load ngsplot/2.47\n";
			print VER "EXEC module load ngsplot/2.47\n";
			print SH "\n# Run default preliminary NGS plot analysis.\n";
			my $ref = $reference{$ip};
			print SH "# \( Full set of ngsplot parameters here: https\:\/\/github.com\/shenlab-sinai\/ngsplot\/wiki\/ProgramArguments101 \)\n";
			print SH "ngs.plot.r -G $ref -FL 150 -L 1000 -R bed -C $outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultConfig.txt -RR 30 -GO none -SC global -VLN 0 -LWD 1 -p 4 -CO blue -O $outputDirectory\/$project\/peaks\/$ip.ngsplot.defaultPlots\n";

			# print SH "\n# Extend peaks from the summit, adding $upstream bp upstream and $downstream bp downstream.\n";
			# print SH "bedtools slop -i $outputDirectory\/$project\/peaks\/$ip.macsPeaks\_summits.bed -g $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes -l $upstream -r $downstream > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.bed\n";
			# print SH "\n# Filter out peaks with an maximum input rpm over 1.\n";
			# print SH "Rscript $NGSbartom/tools/bin/filterOutHighInputPeaks.R  --inputfile=$outputDirectory\/$project\/tracks\/$input.bw --bedfile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.bed --maxInput=1\n";
			# print SH "\n# Take the top peaks (at most 5000) and continue with them.\n";
			# print SH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.bed | head -n 5000 > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.top5k.bed\n";
			# print SH "\n# Make a heatmap and meta plot for expanded peaks.\n";
			# print SH "Rscript $NGSbartom/tools/bin/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.expanded.$upstream.$downstream.max1.top5k.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";

			# print SH "\n# Find Top 100 TSS-proximal peaks, find coordinates of regions around associated TSS's  and make a heatmap and meta plot.\n";
			# print SH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.macsPeaks.anno.txt | awk \'\$10 < $distToTSS {print \$5,\$10,\$11}\' | awk \'\{print \$3\}\' | sort | uniq | head -n 100 > $outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.txt\n";
			# print SH "Rscript $NGSbartom/tools/bin/fromGeneListToTSSbed.R --txdbfile=$txdbfile{$reference{$project}} --assembly=$reference{$project} --geneList=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.txt --up=$upstream --down=$downstream --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			# print SH "Rscript $NGSbartom/tools/bin/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.macsPeaks.100mostOccTSS.$upstream.$downstream.tss.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			if ($buildNGSplot == 1){
			    print SH "module load ngsplot/2.47\n";
			    print VER "EXEC module load ngsplot/2.47\n";
			    print SH "\n# Setting up NGS plot pipeline to run logFC clustered metaPeakPlot.\n";
			    print SH "# For full description of options see https://github.com/ebartom/NGSbartom/blob/master/NGSplotPipeline/NGSplotPipeline.presentation.pdf\n";
			    print SH "# Heatmaps will be calculated for each bed file in the $outputDirectory/$project/analysis/$ip.bedList.txt file ; default is to do each bed file separately.\n";
			    print SH "echo \"$ip\t$outputDirectory\/$project\/peaks\/$ip.macsPeaks.bed\" | cat \>  $outputDirectory/$project/analysis/$ip.bedList.txt\n";
			    print SH "perl $NGSbartom/tools/bin/NGSplotPipeline/makeNGSplots.pl -hss $outputDirectory/$project/scripts/$ip.homer.logFC.metaPeakPlot.clustered.sh -hss2 $outputDirectory/$project/scripts/$ip.homer2.logFC.metaPeakPlot.clustered.sh \\\n";
			    print SH "\t-os $outputDirectory/$project/scripts/$ip.analysis.logFC.metaPeakPlot.clustered.sh \\\n";
			    print SH "\t-o $outputDirectory/$project/NGSplotPipeline/$ip.logFC.metaPeakPlot.clustered \\\n";
			    print SH "\t-vl 0 -lw 1 -g $reference{$project} -fl 150 -p $numProcessors \\\n";
			    print SH "\t-c $ngsFCcomparison \\\n";
			    print SH "\t-csb 0 -ccenb 0 -cepw 1 -cbl 5000 -chs -2,2 -chc blue:white:red -cbo km -cbc 4 -ccd 1 \\\n";
			    print SH "\t-cb $outputDirectory/$project/analysis/$ip.bedList.txt \n";
			    print SH "\n# Run created NGSplot analysis script.\n";
			    print SH ". $outputDirectory/$project/scripts/$ip.analysis.logFC.metaPeakPlot.clustered.sh\n";
			}
			close SH;

		    } elsif ($peakType eq "broad"){
			my (%modulesLoaded, $moduleText);
			$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
			if ($moduleText ne ""){ print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
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
			$modulesLoaded{"bedtools/2.17.0"} = 0;
			$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
			if ($moduleText ne ""){ print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
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
			    print BSH "$NGSbartom/tools/bin/bedToBigBed $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.capped.bed $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bb\n";
			    print BSH "date\n";
			    if ($uploadASHtracks == 1){
				$moduleText = &checkLoad("python/anaconda",\%modulesLoaded);
				if ($moduleText ne ""){ print BSH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"python/anaconda"} = 1;}
				print BSH "\n# Copy bigBed to Amazon S3, for UCSC genome browser to access.\n";
				print BSH "# Note that these files are not visible to browser unless you \"make public\" from within the S3 interface\n";
				print BSH "aws s3 cp $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bb s3://$s3path/$scientist.$project/ --region us-west-2\n";
				if ($modulesLoaded{"python/anaconda"} == 1){ print BSH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
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
			if ($modulesLoaded{"python/anaconda"} == 1){ print BSH "module unload python/anaconda\n"; $modulesLoaded{"python/anaconda"} = 0;}
			print BSH "\nmodule load mpi/openmpi-1.6.3-gcc-4.6.3\nmodule load R/3.2.2\n";
			print BSH "Rscript $NGSbartom/tools/bin/addGenesToBed.R --peakFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed --outputDirectory=$outputDirectory\/$project\/peaks --assembly=$reference{$project} --txdbfile=$txdbfile{$reference{$project}}\n";
			print BSH "date\n";
			if ($buildNGSplot == 1){
			    print BSH "\n# Setting up NGS plot pipeline to run logFC clustered metaPeakPlot.\n";
			    print BSH "# For full description of options see https://github.com/ebartom/NGSbartom/blob/master/NGSplotPipeline/NGSplotPipeline.presentation.pdf\n";
			    print BSH "# Heatmaps will be calculated for each bed file in the $outputDirectory/$project/analysis/$ip.bedList.txt file ; default is to do each bed file separately.\n";
			    print BSH "echo \"$ip\t$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed\n\" | cat \>  $outputDirectory/$project/analysis/$ip.bedList.txt\n";
			    print BSH "perl $NGSbartom/tools/bin/NGSplotPipeline/makeNGSplots.pl -hss $outputDirectory/$project/scripts/$ip.homer.logFC.metaPeakPlot.clustered.sh -hss2 $outputDirectory/$project/scripts/$ip.homer2.logFC.metaPeakPlot.clustered.sh \\\n";
			    print BSH "\t-os $outputDirectory/$project/scripts/$ip.analysis.logFC.metaPeakPlot.clustered.sh \\\n";
			    print BSH "\t-o $outputDirectory/$project/NGSplotPipeline/$ip.logFC.metaPeakPlot.clustered \\\n";
			    print BSH "\t-vl 0 -lw 1 -g $reference{$project} -fl 150 -p $numProcessors \\\n";
			    print BSH "\t-c $ngsFCcomparison \\\n";
			    print BSH "\t-csb 0 -ccenb 0 -cepw 1 -cbl 5000 -chs -2,2 -chc blue:white:red -cbo km -cbc 4 -ccd 1 \\\n";
			    print BSH "\t-cb $outputDirectory/$project/analysis/$ip.bedList.txt \\\n";
			    print BSH "\n# Run created NGSplot analysis script.\n";
			    print BSH ". $outputDirectory/$project/scripts/$ip.analysis.logFC.metaPeakPlot.clustered.sh\n";
			}
			# print BSH "\n# Find summits of peaks.\n";
			# print BSH "Rscript $NGSbartom/tools/bin/fromBedPlusBWtoSummit.R --bedfile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.bed --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			# print BSH "module load bedtools/2.17.0\n";
			# print BSH "bedtools slop -i $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.summits.bed -g $NGSbartom/anno/chromSizes/$reference{$project}\.chrom.sizes -l $upstream -r $downstream > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.bed\n";
			# print BSH "\n# Filter out peaks with an maximum input rpm over 1.\n";
			# print BSH "Rscript $NGSbartom/tools/bin/filterOutHighInputPeaks.R  --inputfile=$outputDirectory\/$project\/tracks\/$input.bw --bedfile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.bed --maxInput=1\n";
			# print BSH "\n# Make a heatmap and meta plot for top 5000 broad peaks.\n";
			# print BSH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.bed | head -n 5000 > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.top5k.bed\n";
			# print BSH "Rscript $NGSbartom/tools/bin/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.expanded.$upstream.$downstream.max1.top5k.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			# print BSH "\n# Find Top 100 TSS-proximal peaks, find coordinates of regions around associated TSS's  and make a heatmap and meta plot.\n";
			# print BSH "sort -nr -k 5 $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.anno.txt | awk \'\$10 < $distToTSS {print \$5,\$10,\$11}\' | awk \'\{print \$3\}\' | sort | uniq | head -n 100 > $outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.txt\n";
			# print BSH "Rscript $NGSbartom/tools/bin/fromGeneListToTSSbed.R --txdbfile=$txdbfile{$reference{$project}} --assembly=$reference{$project} --geneList=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.txt --up=$upstream --down=$downstream --bwfile=$outputDirectory\/$project\/tracks\/$ip.bw\n";
			# print BSH "Rscript $NGSbartom/tools/bin/fromBedPlusBWsToCDTnPlot.R --bedFile=$outputDirectory\/$project\/peaks\/$ip.sicerPeaks.100mostOccTSS.$upstream.$downstream.tss.bed --bwlist=$outputDirectory\/$project\/tracks\/bwlist.txt\n";
			# print BSH "\nmodule unload R\nmodule unload mpi\nmodule load python/anaconda\n";
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
	    my $result = "";
	    if ($scheduler eq "MOAB"){
		$result = `find $outputDirectory/$project/scripts/ -iname \"*callPeaks.sh\" -exec msub {} ./ \\\;`;
		$result =~ s/\s+/\:/g;
		$result =~ s/^://;
		$result =~ s/:$//;
	    } elsif ($scheduler eq "SLURM"){
		$result = `find $outputDirectory/$project/scripts/ -iname \"*callPeaks.sh\" -exec sbatch {} ./ \\\;`;
		$result =~ s/Submitted batch job //g;
		$result =~ s/\s+/\:/g;
		$result =~ s/^://;
		$result =~ s/:$//;
	    }
	    &datePrint( "Need to wait for the following jobs to finish:\n$result");
	    open (SH, ">$outputDirectory/$project/scripts/PeakCallingDependentScript.sh");
	    my (%modulesLoaded, $moduleText);
	    print SH $header;
	    if ($scheduler eq "MOAB"){
		print SH "#MSUB -W depend=afterok:$result\n";
		print SH "#MSUB -N CheckingPeakCallerProgress\n";
		print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
	    } elsif ($scheduler eq "SLURM"){
		print SH "#SBATCH --job-name=CheckingPeakCallerProgress\n";
		print SH "#SBATCH --nodes=1\n";
		print SH "#SBATCH -n $numProcessors\n";
		print SH "#SBATCH --dependency=afterok:$result\n";
	    }
	    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
	    print SH "\necho \"Peaking calling jobs $result have finished.\"\n";
	    close SH;
	    &datePrint("Creating dependent job that will only run after peak callers finish.");
	    my $result2 = "";
	    if ($scheduler eq "MOAB"){
		$result2 = `msub $outputDirectory/$project/scripts/PeakCallingDependentScript.sh`;
		$result2 =~ s/\s+//g;
	    } elsif ($scheduler eq "SLURM"){
		$result2 = `sbatch $outputDirectory/$project/scripts/PeakCallingDependentScript.sh`;
		$result2 =~ s/Submitted batch job //g;
		$result2 =~ s/\s+//g;
	    }
	    my $jobfinished = "no";
	    # Wait until the job is Complete.
	    &datePrint("Waiting for job $result2 to finish. (each . = 300 seconds)");
	    waitForJob($result2);
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
		my (%modulesLoaded, $moduleText);
		print SH $header;
		if ($scheduler eq "MOAB"){
		    print SH "#MSUB -N $peakset\_diffPeak\n";
		    print SH "#MSUB -l nodes=1:ppn=$numProcessors\n";
		} elsif ($scheduler eq "SLURM"){
		    print SH "#SBATCH --job-name=$peakset\_diffPeak\n";
		    print SH "#SBATCH --nodes=1\n";
		    print SH "#SBATCH -n $numProcessors\n";
		}
		print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
		$moduleText = &checkLoad("bedtools/2.17.0",\%modulesLoaded);
		if ($moduleText ne ""){ print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"bedtools/2.17.0"} = 1;}
		$moduleText = &checkLoad("samtools/1.6",\%modulesLoaded);
		if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.6"} = 1;}
		#$moduleText = &checkLoad("samtools/1.2",\%modulesLoaded);
		#if ($moduleText ne ""){	print SH $moduleText; print VER "EXEC $moduleText"; $modulesLoaded{"samtools/1.2"} = 1;}
		print SH "module load R/3.2.2\n";
		print VER "EXEC module load R/3.2.2\n";
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
		    print SH "cat @samples | \\\nawk \'\{printf \"\%s\\t\%d\\t\%d\\n\",\$1,\$2,\$3\}\' - | \\\n";
		    print SH "\tbedtools sort -i - | \n bedtools merge -i - | \\\n\t";
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
		system("sh $NGSbartom/tools/bin/moduleLoadSamtools.sh");
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
		print SH "perl $NGSbartom\/tools\/bin\/makePeakCountsTable.pl $outputDirectory\/$project\/peaks\/ $outputDirectory\/$project\/analysis\/ $peakset\n";
		print SH "\n# Generate MDS plot for peakset $peakset\n";
		print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --outputDirectory=$outputDirectory\/$project\/analysis/ --countFile=$outputDirectory\/$project\/analysis\/$peakset.all.counts.txt --numCores=$numProcessors --runMDS=1\n";
		&datePrint("Looking for $outputDirectory\/$project\/$peakset.comparisons.csv");
		if (-e "$outputDirectory\/$project\/$peakset.comparisons.csv"){
		    print SH "\n# Find Differential Peaks for peakset $peakset\n";
		    print SH "Rscript $NGSbartom/tools/bin/runEdgeRrnaSeq.2.R --assembly=$reference{$project} --outputDirectory=$outputDirectory\/$project\/analysis/ --countFile=$outputDirectory\/$project\/analysis\/$peakset.all.counts.txt --numCores=$numProcessors --comparisonFile=$outputDirectory\/$project\/$peakset.comparisons.csv --runMDS=0\n";

		    ## For each comparison, go through and output Bed file of significant peaks
		} else {
		    &datePrint("To find differentially expressed peaks in peakset $peakset, create a comparisons file at $outputDirectory\/$project\/$peakset.comparisons.csv")
		}
		print SH "date\n";
		close(SH);
		if ($runDiffPeaks == 1){
		    &datePrint("Starting job $outputDirectory/$project/scripts/$peakset\_diffPeaks.sh");
		    if ($scheduler eq "MOAB"){
			`msub $outputDirectory/$project/scripts/$peakset\_diffPeaks.sh`;
		    }elsif ($scheduler eq "SLURM"){
			`sbatch $outputDirectory/$project/scripts/$peakset\_diffPeaks.sh`;
		    }
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
    #    print SH "\n#If there are any modules loaded, remove them.\nmodule purge\n\n";
    # 	    print SH "\necho \"Peaking calling jobs $result have finished.\"\n";
    # 	    close SH;
    # 	    &datePrint("Creating dependent job that will only run after diff peak scripts finish.");
    # 	    my $result2 = `msub $outputDirectory/$project/scripts/diffPeakDependentScript.sh`;
    # 	    $result2 =~ s/\s+//g;
    # 	    my $jobfinished = "no";
    # 	    # Wait until the job is Complete.
    # 	    &datePrint("Waiting for job $result2 to finish. (each . = 300 seconds)");
    # 	    # Check qstat every 300 seconds, adding a "." every time you check.
    # 	    until ($jobfinished eq "Completed"){
    # 		sleep(300);
    # 		$jobfinished = `checkjob $result2 | grep ^State:`;
    # 		print STDERR ".";
    # 		if ($jobfinished =~ /State: (\w+)\s+/){ $jobfinished = $1;}
    # 		print STDERR "$jobfinished";
    # 	    }
    # 	    print STDERR "\n";
    # 	    &datePrint("Job $result2 done.  Continuing.");
    #}
    #}
}


close(VER);
if (-e "$outputDirectory\/metadata\/Ceto.run.$type.$timestamp.txt"){
    `sort $outputDirectory\/metadata\/Ceto.run.$type.$timestamp.txt | uniq > $outputDirectory\/metadata\/Ceto.run.$type.$timestamp.uniq.txt`;
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

sub waitForJob {
	my $jobId = $_[0];
	my $qstat_output = '';
	if ($scheduler eq "MOAB"){
	    until ($qstat_output =~ /Unknown Job Id/ || $qstat_output =~ /job_state = C/) {
		sleep(300);
		$qstat_output = `qstat -f $jobId.qsched03.quest.it.northwestern.edu 2>&1`;
		print STDERR ".";
	    }
	} elsif ($scheduler eq "SLURM"){
	    until ($qstat_output =~ /COMPLETED/){
		sleep(300);
		$qstat_output = `sacct -j $jobId -n -b`;
		print STDERR ".";
#		print STDERR "$qstat_output";
	    }
	}
	print STDERR "Job $jobId is complete, scheduler $scheduler\n";
	if ($scheduler eq "SLURM"){
	    `seff $jobId > $outputDirectory/metadata/SLURM.checkJob.$jobId.seff.txt`;
	}
}

sub checkLoad {
    my $module = $_[0];
    my $hashref = $_[1];
    my %modulehash = %$hashref;
#    print STDERR "# Check load!\n";
#    print STDERR "$module\n";
#    print STDERR "$modulehash{$module}\n";
    if ($modulehash{$module} == 1){
	# Skip these module loads.
#	print STDERR "skip load\n";
	return ("");
    } elsif ((!exists($modulehash{$module})) || ($modulehash{$module} == 0)){
#	print STDERR "must load\n";
	return ("module load $module\n");
    }
}
