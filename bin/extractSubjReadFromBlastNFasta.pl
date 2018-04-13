#!/usr/local/bin/perl -w

use strict;

if($#ARGV < 0) {
  print $#ARGV;
  print STDERR "Useage: $0 < Blast file > < Fasta file > [ min percent identity ]\n";
  exit;
}

my $thresh = 0;
my $minLength = 16;
my $blastfile = $ARGV[0];
my $fastafile = $ARGV[1];
$thresh = $ARGV[2];
my %blasthit;

open(BL,$blastfile);
while (<BL>){
    chomp $_;
    # Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    if ($_ !~ /^\#/){
	my ($query,$subj,$perID,$alnLength,$mismatch,$gapOpens,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split (/\t/,$_);
	#print STDERR "$subj\t$query\n";
	if ($perID >= $thresh){
	    $blasthit{$subj} = $_;
	}
    }
}

my $string = "";
my $header = "";
open(FA,$fastafile);
while(<FA>){
    chomp $_;
    $_ =~ s/\r\n/\n/g;
    $_ =~ s/\r//g;
    if ($_ =~ /^([ACTGN]+)$/){
	chomp $_;
	$string .= $1;
#	print STDERR "Found seq $1\t$header\n";
    } elsif ($_ =~ /^>([\.\w\:\=\_\-]+)/){
	if ($string ne ""){
	    if (exists($blasthit{$header})){
		my ($query,$subj,$perID,$alnLength,$mismatch,$gapOpens,$qstart,$qend,$sstart,$send,$evalue,$bitscore) = split (/\t/,$blasthit{$header});
		my $chimera = substr($string,$send);
		if (length($chimera) > $minLength){
		    print ">$query.$subj.chimera\n$chimera\n";
		}
	    }
#	    else { print STDERR "$header\n";}
	    $string = "";
	}
	$header = $1;
    }
}

