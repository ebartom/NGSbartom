#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < Seqfile > [ gene list ]\n";
  exit;
}

#
# This script takes a fasta formatted file (where the sequences
# can be interupted by carriage returns) and identifies all sequences
# with the same gene name.  It returns all of the sequences in fasta 
# format (with no carriage returns).  Optionally it also takes a gene 
# list, and returns only those genes.
#
# This was intially developed to pull out the longest 3' utr sequence
# for each human gene.
#

my $seqfile = $ARGV[0];
my $genelist = $ARGV[1];

open (FA,$seqfile);
my $seq = "";
my $header = "";
my %longestUTR;
my %lengths;
my %seqs;
my %cgns;
my %transcripts;
while (<FA>){
  if ($_ =~ /^>([\w\s\d\-\|\_\.\'\:\;]+)\n/){
 #   print STDERR ("Found new sequence $1\n");
    if (($header ne "") || ($seq eq "SEQUENCEUNAVAILABLE")){
  #    print STDERR ("Processing old sequence $header\n");
      $seq = uc($seq);
      $seq =~ s/\s//g;
      my ($txID,$cgn,@stuff) = split(/\|/,$header);
      $cgns{$txID} = $cgn;
      $seqs{$txID} = $seq;
      if (exists($transcripts{$cgn})){
	  $transcripts{$cgn} .= ",$txID";
      } else {
	  $transcripts{$cgn} = $txID;
      }
    }
    $header = $1;
    $seq = "";
  } elsif ($_ !~ /Sequence unavailable/) {
    $seq .= $_;
  }
}

my %whitelist;
if ($genelist){
    open (IN,$genelist);
    while (<IN>){
	chomp $_;
	$whitelist{$_} = "1";
    }
}


my ($txID,$cgn,@stuff);
if (($header ne "") || ($seq eq "SEQUENCEUNAVAILABLE")){
#  print STDERR ("Processing old sequence $header\n");
  $seq = uc($seq);
  $seq =~ s/\s//g;
  ($txID,$cgn,@stuff) = split(/\|/,$header);
  $cgns{$txID} = $cgn;
  $seqs{$txID} = $seq;
  if (exists($transcripts{$cgn})){
      $transcripts{$cgn} .= ",$txID";
  } else {
      $transcripts{$cgn} = $txID;
  }
}

my @cgn = keys(%transcripts);

foreach my $cgn (@cgn){
    my @transcripts = split (/,/,$transcripts{$cgn});
    foreach my $txID (@transcripts){
	my $seq = $seqs{$txID};
	if (length($seq) > 0){
	    if ($genelist){
		if (exists($whitelist{$cgn})){
		    print ">$txID\|$cgn\n$seq\n";
		}
	    }else {
		print ">$txID\|$cgn\n$seq\n";
	    }
	}
    }
}
