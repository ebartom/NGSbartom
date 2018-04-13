#!/usr/bin/perl -w
use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $samfile = $ARGV[0];

if ($#ARGV <0){
    die "Usage:  $0 <samfile>\n";
}

my (%kmers);
open (SAM,$samfile);
while(<SAM>){
    chomp $_;
    my ($name,$flag,$chr,$start,$mapq,$cigar,$star,$zero,$zero2,$seq,$qual,@other) = split(/\t/,$_);
    my $string = "$chr:$start:$cigar";
    if ($string ne "*:0:*"){
	if (exists($kmers{$seq})){
	    $kmers{$seq} .= "|$string";
	} else {
	    $kmers{$seq} = $string;
	}
    }
}

my $bedfile = $samfile;
$bedfile =~ s/.sam$/.bed/g;

open(BED,">$bedfile");

foreach my $kmer (keys(%kmers)){
    my $string = $kmers{$kmer};
#    print STDERR "$string\n";
    my @matches = split(/\|/,$string);
    my @uniqMatches = uniq(@matches);
    $string = join("\|", @uniqMatches);
    my $numHits = $#uniqMatches + 1;
    print "$kmer\t$numHits\t$string\n";
    my $k = length($kmer);
    for (my $i=0; $i < $numHits;$i++){
	my $match = $uniqMatches[$i];
	my ($chr,$start,$cigar) = split(/\:/,$match);
	my $stop = $start + $k;
	print BED "$chr\t$start\t$stop\t$kmer.m$i\n";
    }
}
