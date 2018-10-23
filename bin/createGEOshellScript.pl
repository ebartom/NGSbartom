#!/usr/bin/perl -w
use List::MoreUtils qw(uniq);

use strict;
use warnings;

my $sampleFile = $ARGV[0];

#-----------------------------------------------------------
# This script takes a tab delimited file in the form of
# OldSampleName\tNewSampleName\tTANGO-ID\tExperimentType\n
#
# It looks for the old samples in the tango directories, and
# copies them into subdirectories in the working directory.
# It concatenates sequencing lanes to create new fastq files
# with all of the reads for a given sample, renaming it with
# the new sample name.  It will also copy over track files
# which can be used as "processed files" in a GEO repository.
#
# All of the necessary commands are written to STDOUT.  I send
# them to a shell script, which I can then run (and edit if
# absolutely necessary).  Sample lines which are preceded by
# a pound sign (#) will be skipped by this perl script.  
#
# Elizabeth Bartom 10/23/2018
# ebartom@northwestern.edu
#------------------------------------------------------------

if ($#ARGV <0){
    die "Usage:  $0 < sampleFile >\n";
}
my %tango;
my %nameChange;
my %expType;

my $path1 =  "/projects/b1025/tango/";
my $path2 =  "/projects/b1025/xvault/tango/";

open(IN,$sampleFile);
while(<IN>){
    chomp $_;
    if (($_ !~ /^#/) && ($_ ne "")) {
	my ($oldName,$newName,$tangoID,$type) = split(/\t/,$_);
	$oldName =~ s/\-/\_/g;
	$nameChange{$oldName} = $newName;
	$tango{$oldName} = $tangoID;
	$expType{$oldName} = uc($type);
    }
}

my @oldSamples = keys(%tango);
my @tangos = uniq(values(%tango));
#my @categories = uniq(values(%dir));

my %tangoShell;
foreach my $sample (@oldSamples){
    my $tango = $tango{$sample};
    my $oldSample = $sample;
    my $type = $expType{$sample};
    $oldSample =~ s/\_/\-/g;
    if (!exists($tangoShell{$tango})){ $tangoShell{$tango} = "";}
    $tangoShell{$tango} .= "\# $sample\t$nameChange{$sample}\t$tango\n";
    my $fastqs = "";
    my $tracks = "";
    my $trackDir = "tracks";
    if ($type =~ /RNA/) { $trackDir = "bam";}
    $fastqs = `ls $path1/$tango/$tango/fastq/$sample*fastq.gz`;
    if ($fastqs eq ""){
	$fastqs = `ls $path1/*$tango*/$tango/fastq/$sample*fastq.gz`;
    } else { $tracks = `ls $path1/$tango/$tango/$trackDir/$oldSample*bw`;}
    if ($fastqs eq ""){
	$fastqs = `ls $path2/$tango/$tango/fastq/$sample*fastq.gz`;
    } elsif ($tracks eq "") { $tracks = `ls $path1/*$tango*/$tango/$trackDir/$oldSample*bw`;}
    if ($fastqs eq ""){
	$fastqs = `ls $path2/*$tango*/$tango/fastq/$sample*fastq.gz`;
    } elsif ($tracks eq "") { $tracks = `ls $path2/$tango/$tango/$trackDir/$oldSample*bw`;}
    if ($fastqs eq ""){
	my $myTango = $tango;
	$myTango =~ s/\-/\-\*/g;
	$fastqs = `ls $path1/*$myTango*/$tango/fastq/$sample*fastq.gz`;
	if ($fastqs eq ""){
	    $fastqs = `ls $path2/*$myTango*/$tango/fastq/$sample*fastq.gz`;
	} elsif ($tracks eq "") { $tracks = `ls $path1/*$myTango*/$tango/$trackDir/$oldSample*bw`;}
	if ($tracks eq "") { $tracks = `ls $path2/*$myTango*/$tango/$trackDir/$oldSample*bw`;}
    } elsif ($tracks eq "") { $tracks = `ls $path2/*$tango*/$tango/$trackDir/$oldSample*bw`;}
    if ($fastqs eq ""){
	$tangoShell{$tango} .= "\# FAILED: $path1/*$tango*/$tango/fastq/$sample*fastq.gz\t$path2/*$tango*/$tango/fastq/$sample*fastq.gz\n";
    }
    $fastqs =~ s/\n/ /g;
    $tracks =~ s/\n/ /g;
    # First Fastqs.
    $tangoShell{$tango} .= "cp $fastqs $tango/originalFastq/\n";
    $tangoShell{$tango} .=  "zcat $tango/originalFastq/$sample*fastq.gz | gzip -c > $tango/renamedFastq/$nameChange{$sample}.fastq.gz\n";
    $tangoShell{$tango} .= "md5sum $tango/renamedFastq/$nameChange{$sample}.fastq.gz\n\n";
    # Then Tracks.
    $tangoShell{$tango} .= "# Copying tracks for experiment type $type.\n";
    $tangoShell{$tango} .= "cp $tracks $tango/originalTracks/\n";
    if ($type =~ /RNA/){
	$tangoShell{$tango} .=  "cp $tango/originalTracks/$oldSample.plus.bw $tango/renamedTracks/$nameChange{$sample}.plus.bw\n";
	$tangoShell{$tango} .=  "cp $tango/originalTracks/$oldSample.minus.bw $tango/renamedTracks/$nameChange{$sample}.minus.bw\n";
	$tangoShell{$tango} .= "md5sum $tango/renamedTracks/$nameChange{$sample}.plus.bw\n";
	$tangoShell{$tango} .= "md5sum $tango/renamedTracks/$nameChange{$sample}.minus.bw\n";
    }else {
	$tangoShell{$tango} .=  "cp $tango/originalTracks/$oldSample.bw $tango/renamedTracks/$nameChange{$sample}.bw\n";
	$tangoShell{$tango} .= "md5sum $tango/renamedTracks/$nameChange{$sample}.bw\n\n";
    }
}

foreach my $tango (@tangos) {
    print "\# Setup $tango\n";
    print "mkdir $tango\n";
    print "mkdir $tango/originalFastq\n";
    print "mkdir $tango/renamedFastq\n";
    print "mkdir $tango/originalTracks\n";
    print "mkdir $tango/renamedTracks\n\n";
    print "\# Copy over fastqs and tracks (raw and processed files)\n";
    print $tangoShell{$tango};
}
