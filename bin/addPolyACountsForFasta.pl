#!/opt/local/bin/perl -w

################################################################
# This script will take an input sequence file in Fasta Format.
# For each sequence in the file, the script will find the longest
# stretch of polyA at least 6 bp in length.  It will return the
# sequence name, the sequence, the sequence length, a flag for
# whether the sequence has a polyA at least 6 bp in length, the
# length of the longest stretch, and the position of the longest
# polyA relative to both the head and tail of the sequence.
#
# Elizabeth Bartom, Northwestern University, December 2017.
#
################################################################

use strict;

if($#ARGV < 0) {
  print $#ARGV;
  print STDERR "Useage: $0 < inputFasta >\n";
  exit;
}

my $k = 6;
my $fastafile = $ARGV[0];

my @alphabet = ("A");

my $string = "";
my $header = "";
my %intro;
my (%Aflag);
my (%Adata);

# Read in each sequence in the sequence file, removing extra whitespace.
# Find and record the longest polyA tract, and whether it is over 6
# nucleotides in length.
open(IN,$fastafile);
while(<IN>){
    chomp $_;
    $_ =~ s/\r\n/\n/g;
    $_ =~ s/\r//g;
    if ($_ =~ /^([ACTGN]+)$/){
	#	print STDERR "Sequence Line $string\n";
	chomp $_;
	$string .= $1;
    } elsif ($_ =~ /^>(.+)$/){
	if ($string ne ""){
	    my $strlen = length($string);
	    $intro{$header} = "$header\t$string\t$strlen";
	    my ($cgn,$txID) = split(/\s+/,$header);
	    if (!exists($Aflag{$cgn})){ 
		$Aflag{$cgn} = 0;
	    }
	    # polyA
	    if ($string =~ s/([CGTN]|\b)(A{$k,})([CGTN]|\b)/${1}${2}${3}/g){
		my $start = $-[2];
		my $end = $+[2];
		$Adata{$header} = "\t$2\t".length($2)."\t".$start."\t";
		$Adata{$header} .= $strlen-$end;
		$Aflag{$cgn} = 1;
	    } else {
		$Adata{$header} = "\t\t0\tna\tna";
	    }
	}
	$header = $1;
	$header =~ s/\|/\t/g;
	$string = "";
    }
}

# Set up the first line of the file with column labels.
print "EnsID\tGene\tSeq\tLength";
foreach my $nuc (@alphabet){
    print "\t".$nuc."Flag\t$nuc\tmaxLength\t5pDist\t3pDist";
}
print "\n";

# Print the data for each sequence.
foreach my $header (keys(%intro)){
    print $intro{$header}."\t";
    my ($cgn,$Txid) = split(/\s+/,$header);
    print $Aflag{$cgn}.$Adata{$header}."\n";
}
