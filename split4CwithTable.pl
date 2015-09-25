#!/usr/bin/perl -w
use strict;
use warnings;
use List::Util qw( min max );

my $fastQdir = $ARGV[0];
my $sampleDescFile = $ARGV[1];
# GACGTAATTTTGGT
my $maxMismatch = $ARGV[2];
my $cutSiteLength = 6;
#my $resSite = "AAGCTT";

if ($#ARGV <1){
    die "Usage:  $0 <fastQdirectory> <sampleDescriptionTable> <maxMismatch (viewpoint)>\n";
}

my %index;
my %viewpoint;
my %sample;
my %reads;
my %adapters;

open(IN,$sampleDescFile);
while(<IN>){
    chomp $_;
    if ($_ !~ /^Sample/){
	print STDERR "$_\n";
	my ($sampleName,$indexSeq,$viewpointSeq,$assembly,$viewpoint) = split(/\t/,$_);
	$adapters{$viewpointSeq} = 1;
	$viewpoint{$sampleName} = $viewpointSeq;
	$index{$sampleName} = $indexSeq;
#	print STDERR "$sampleName\t$indexSeq,$viewpointSeq\n";
	$sample{$indexSeq}{$viewpointSeq} = $sampleName;
	$reads{$indexSeq}{$viewpointSeq} = "";
    }
}

my $fastqFiles = `ls $fastQdir\/Undetermined*.fastq`;
my @fastqFiles = split(/\s+/,$fastqFiles);
print STDERR "FastQfiles:  @fastqFiles\n";

foreach my $seqfile (@fastqFiles){
#    open (OUT,">$seqfile.$adapter.catch.m$maxMismatch.fq");
    
    foreach my $adapter (keys(%adapters)){
	open (FQ,$seqfile);
	print STDERR "Sequence file: $seqfile\n";
	my $resSite = substr($adapter,0-$cutSiteLength,$cutSiteLength);
	my $longAdapter = $adapter;
	$adapter = substr($adapter,0,length($adapter)-$cutSiteLength);
	print STDERR "Restriction enzyme cut site = \"$resSite\"\n";
	print STDERR "Viewpoint adapter sequence = \"$adapter\"\n";
	my $mis_adapter = 0;
	my $len_adapter = length($adapter);
	my $lineCount = 5;
	my ($header,$seq,$qual,$header2);
	while(<FQ>){
	    chomp $_;
	    #  print STDERR "$lineCount\t$_\n";
	    if (($_ =~ /^\@/)&& ($lineCount == 5)){ $header =$_; $lineCount = 1;}
	    elsif (($lineCount == 2)){ $seq = $_;}
	    elsif (($lineCount == 3)){ $header2 = $_}
	    elsif ($lineCount == 4){ 
	   $qual = $_;
	   #    print STDERR "$header $seq\n";	       
	   my $primer = "";
	   my $catch = "";
	   my $type = "??";
	   my $catchStart = 0;
	   my $catchQual = "";
	   my $index = "";
	   if ($seq =~ /([\w\.]+)$resSite([\w\.]+)/){
	       $catchStart = $-[2];
	       $catchQual = substr($qual,$catchStart-length($resSite));
	       $primer = $1;
	       $catch = $2;
#	            print STDERR "A: $index $primer $resSite $catch\n";
	       my $primerLen = length($primer);
	       my $catchLen = length($catch);
	       my $numN = () = $catch =~ /N/gi;
	       #      print STDERR "$primer  $primerLen\t$catch $catchQual\n";
	       if (($primerLen >= $len_adapter)&& ($catchLen >= 10) && ($numN < ($catchLen/2))){
		   my $A = substr($primer,length($primer)-$len_adapter,$len_adapter);
		   my $B = $adapter;
		   $index = substr($primer,0,(length($primer)-$len_adapter));
#		   print STDERR "$seq\n";
#		   print STDERR "B: \"$index\" \"$adapter\" $A $resSite $catch\n";
		   $mis_adapter = &mismatch_count($A,$B);
#		   print STDERR "\"$A\" vs. \"$B\"\t$mis_adapter\n";
		   if ($mis_adapter <= $maxMismatch) {
		       if (exists($sample{$index}{$longAdapter})){
			   my $sampleName = $sample{$index}{$longAdapter};
#			   print STDERR "$sampleName\t$index\t$adapter\t$longAdapter\n";
			   $reads{$index}{$longAdapter} .= "$header $sampleName $mis_adapter $numN\n$resSite$catch\n+\n$catchQual\n";
			   #  $type = "Adapter"; print OUT "$header $index $mis_adapter $numN\n$resSite$catch\n+\n$catchQual\n";
		       }
		   }
	       } 
	   }
	    }
	    #print STDERR "$lineCount\n";
	    $lineCount++;
	}
	close(FQ);
    }
}

foreach my $sample (keys(%index)){
    my $outfile = "$fastQdir\/$sample.m$maxMismatch.fastq";
    open(OUT,">$outfile");
    my $index = $index{$sample};
    my $adapter = $viewpoint{$sample};
    print STDERR "Printing reads for $sample to $outfile\n";
    print OUT $reads{$index}{$adapter};
    close OUT;
    my $cmd = "gzip $fastQdir\/$sample.m$maxMismatch.fastq";
    system($cmd);
}

sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }
    
