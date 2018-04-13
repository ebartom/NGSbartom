#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print $#ARGV;
  print STDERR "Useage: $0 < kmer length > < inputFasta1 > < inputFasta2 > < pseudocount >\n";
  exit;
}

my $k = $ARGV[0];
my $fastafile1 = $ARGV[1];
my $fastafile2 = $ARGV[2];
my $pseudocount = $ARGV[3];

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

my %kmerspaces1;
my %kmerspaces2;

foreach my $kmer (@kmerlist){
  if ($kmer ne ""){
    if (!exists($seen{$kmer})){
      push(@cleanList,$kmer);
      $seen{$kmer} = 1;
      $kmerspaces1{$kmer} = 0;
      $kmerspaces2{$kmer} = 0;
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
my %freqs1;
my %freqs2;
my $gene1num = 0;
print STDERR "Starting first file: $fastafile1\n";
open(IN1,$fastafile1);
while(<IN1>){
    chomp $_;
    $_ =~ s/\r\n/\n/g;
    $_ =~ s/\r//g;
    if ($_ =~ /^([ACTGN]+)$/){
	chomp $_;
	$string .= $1;
    } elsif ($_ =~ /^>(.+)$/){
	if (($string ne "") && length($string) >= $k){
	    $gene1num ++;
	    #	    print "$header\t$string\t".length($string);
	    foreach my $kmer (@kmerlist){
		#		my $kmercount = () = $string =~ /$kmer/g;
		my $kmercount = () = $string =~ /(?=\Q$kmer\E)/g;	
		my $kmerspace = length($string)-$k+1;
		$freqs1{$kmer} += ($kmercount+$pseudocount)/($kmerspace+$pseudocount);
	    }
	    #	    print "\n";
	}
	$header = $1;
	$header =~ s/\|/\t/g;
	$string = "";
    }
}
print STDERR "Starting second file: $fastafile2\n";
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
	if (($string ne "") && length($string) >= $k){
	    $gene2num ++;
	    #	    print "$header\t$string\t".length($string);
	    foreach my $kmer (@kmerlist){
		#		my $kmercount = () = $string =~ /$kmer/g;
		my $kmercount = () = $string =~ /(?=\Q$kmer\E)/g;	
		my $kmerspace = length($string)-$k+1;
		$freqs2{$kmer} += ($kmercount+$pseudocount)/($kmerspace+$pseudocount);
	    }
	}
	$header = $1;
	$header =~ s/\|/\t/g;
	$string = "";
    }
}

print "Kmer\tAverage KmerFreqs $fastafile1 $gene1num\tAverage KmerFreqs $fastafile2 $gene2num\tTI ($fastafile1 KmerFreqs / Num Seqs)/($fastafile2 KmerFreqs / Num Seqs)\n";
foreach my $kmer (@kmerlist){
    my $kmerfreq1 = ($freqs1{$kmer})/$gene1num;
    my $kmerfreq2 = ($freqs2{$kmer})/$gene2num;
    my $TI = 0;
    if ($kmerfreq2 > 0){
	$TI = $kmerfreq1/$kmerfreq2;
    }
    print "$kmer\t$kmerfreq1\t$kmerfreq2\t$TI\n";
}
