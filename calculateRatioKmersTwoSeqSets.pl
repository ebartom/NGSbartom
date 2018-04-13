#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print $#ARGV;
  print STDERR "Useage: $0 < kmer length > < inputFasta1 > < inputFasta2 >\n";
  exit;
}

my $k = $ARGV[0];
my $fastafile1 = $ARGV[1];
my $fastafile2 = $ARGV[2];

my @alphabet = ("A","C","G","T");
#print "$k @alphabet\n";

my @kmerlist = @alphabet;
for (my $i = 1; $i < $k; $i++) {
#  print "$i\n";
  my @newlist;
  foreach my $word (@kmerlist){
    foreach my $letter (@alphabet){
      my $newWord = $word.$letter;
      #      print "$word + $letter = $newWord\n";
      push (@newlist, $newWord);
    }
  }
  @kmerlist = @newlist;
#  push(@kmerlist,@newlist);
}
print STDERR "Kmers @kmerlist\n";
my @cleanList;
my %seen;

foreach my $kmer (@kmerlist){
  if ($kmer ne ""){
    if (!exists($seen{$kmer})){
      push(@cleanList,$kmer);
      $seen{$kmer} = 1;
    }
  }
}
@kmerlist = @cleanList;
print STDERR "KmersClean @kmerlist\n";
my $numKmers = $#kmerlist +1;
print STDERR "Found $numKmers Kmers.\n";

my $string = "";
my $header = "";
#print "EnsID\tGene\tSeq\tLength";
#foreach my $kmer (@kmerlist){
#    print "\t$kmer";
#}
#print "\n";
my %counts1;
my %counts2;
my $gene1num = 0;
open(IN1,$fastafile1);
while(<IN1>){
    chomp $_;
    $_ =~ s/\r\n/\n/g;
    $_ =~ s/\r//g;
    if ($_ =~ /^([ACTGN]+)$/){
	chomp $_;
	$string .= $1;
    } elsif ($_ =~ /^>(.+)$/){
	if ($string ne ""){
	    $gene1num ++;
#	    print "$header\t$string\t".length($string);
	    foreach my $kmer (@kmerlist){
#		my $kmercount = () = $string =~ /$kmer/g;
		my $kmercount = () = $string =~ /(?=\Q$kmer\E)/g;
		$counts1{$kmer} += $kmercount;
	    }
#	    print "\n";
	}
	$header = $1;
	$header =~ s/\|/\t/g;
	$string = "";
    }
}
open(IN2,$fastafile2);
my $gene2num = 0;
while(<IN2>){
    chomp $_;
    $_ =~ s/\r\n/\n/g;
    $_ =~ s/\r//g;
    if ($_ =~ /^([ACTGN]+)$/){
	chomp $_;
	$string .= $1;
    } elsif ($_ =~ /^>(.+)$/){
	if ($string ne ""){
	    $gene2num ++;
#	    print "$header\t$string\t".length($string);
	    foreach my $kmer (@kmerlist){
#		my $kmercount = () = $string =~ /$kmer/g;
		my $kmercount = () = $string =~ /(?=\Q$kmer\E)/g;
		$counts2{$kmer} += $kmercount;
	    }
#	    print "\n";
	}
	$header = $1;
	$header =~ s/\|/\t/g;
	$string = "";
    }
}

print "Kmer\t$fastafile1 KmerCounts\t$fastafile1 SeqNum\tNormCounts\t$fastafile2 KmerCounts\t$fastafile2 SeqNum\tNormCounts\tTI ($fastafile1 NormCounts)/($fastafile2 NormCounts)\n";
foreach my $kmer (@kmerlist){
    my $ratio1 = $counts1{$kmer}/$gene1num;
    my $ratio2 = $counts2{$kmer}/$gene2num;
    my $TI = 0;
    if ($ratio2 > 0){
	$TI = $ratio1/$ratio2;
    }
    print "$kmer\t$counts1{$kmer}\t$gene1num\t$ratio1\t$counts2{$kmer}\t$gene2num\t$ratio2\t$TI\n";
}
