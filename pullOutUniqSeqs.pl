#!/opt/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print STDERR "Useage: $0 < Seqfile >\n";
  exit;
}

#
# This script takes a fasta formatted file (where the sequences
# can be interupted by carriage returns) and identifies all sequences
# with the same sequence.  It returns a unique set, with the count
# added to the header.
#
# This was developed in order to pull out a unique set of read sequences.
#

my $seqfile = $ARGV[0];

open (FA,$seqfile);
my $seq = "";
my $header = "";
my %seqs;
my %headers;
while (<FA>){
  if ($_ =~ /^>([\w\s\d\-\=\|\_\.\'\:\;]+)\n/){
#    print STDERR ("Found new sequence $1\n");
    if (($header ne "") || ($seq eq "SEQUENCEUNAVAILABLE")){
#      print STDERR ("Processing old sequence $header\n");
      $seq = uc($seq);
      $seq =~ s/\s//g;
      $seqs{$header} = $seq;
      if (!exists($headers{$seq})){
	  $headers{$seq} = $header;
      } else {
	  $headers{$seq} .= "|$header";
      }
    }
    $header = $1;
    $seq = "";
  } elsif ($_ !~ /Sequence unavailable/) {
    $seq .= $_;
  }
}

foreach my $seq (keys(%headers)){
    my @headers = split(/\|/,$headers{$seq});
    my $count = $#headers + 1;
    print ">$headers[0] $count\n$seq\n";
}
