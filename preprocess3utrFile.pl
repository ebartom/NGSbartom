#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < Seqfile >\n";
  exit;
}

#
# This script takes a fastq formatted file (where the sequences
# can be interupted by carriage returns) and identifies all sequences
# with the same gene name.  It returns only the longest gene name, 
# sequence pair, in fasta format (with no carriage returns).
#
# This was intially developed to pull out the longest 3' utr sequence
# for each human gene.
#

my $seqfile = $ARGV[0];

open (FA,$seqfile);
my $seq = "";
my $header = "";
my %longestUTR;
my %lengths;
my %seqs;
while (<FA>){
  if ($_ =~ /^>([\w\s\d\-\|\_\.\'\:\;]+)\n/){
    print STDERR ("Found new sequence $1\n");
    if (($header ne "") || ($seq eq "SEQUENCEUNAVAILABLE")){
      print STDERR ("Processing old sequence $header\n");
      $seq = uc($seq);
      $seq =~ s/\s//g;
      my ($txID,$cgn,@stuff) = split(/\|/,$header);
      $lengths{$txID} = length($seq);
      $seqs{$txID} = $seq;
      if (!exists($longestUTR{$cgn})){
	$longestUTR{$cgn}= $txID;
      } elsif ($lengths{$longestUTR{$cgn}} < length($seq)){
	$longestUTR{$cgn}= $txID;
      }
    }
    $header = $1;
    $seq = "";
  } elsif ($_ !~ /Sequence unavailable/) {
    $seq .= $_;
  }
}

my ($txID,$cgn,@stuff);
if (($header ne "") || ($seq eq "SEQUENCEUNAVAILABLE")){
  print STDERR ("Processing old sequence $header\n");
  $seq = uc($seq);
  $seq =~ s/\s//g;
  ($txID,$cgn,@stuff) = split(/\|/,$header);
  $lengths{$txID} = length($seq);
  $seqs{$txID} = $seq;
  if (!exists($longestUTR{$cgn})){
    $longestUTR{$cgn}= $txID;
  } elsif ($lengths{$longestUTR{$cgn}} < length($seq)){
    $longestUTR{$cgn}= $txID;
  }
}

my @cgn = keys(%longestUTR);

foreach my $cgn (@cgn){
  my $txID = $longestUTR{$cgn};
  my $seq = $seqs{$txID};
  if (length($seq) > 0){
    print ">$txID\|$cgn\n$seq\n";
  }
}
