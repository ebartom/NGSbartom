#!/usr/bin/perl -w
use strict;
use warnings;

my $fastQdir = $ARGV[0];
my $outputDir = $ARGV[1];
my $projectID = $ARGV[2];
my $sampleSheet = $ARGV[3];

if ($#ARGV <0){
    die "Usage:  $0 <fastQdirectory> <outputDirectory> <projectID|MOLNG> <sample_report>\n";
}

if ($sampleSheet eq ""){
    $sampleSheet = "$fastQdir\/Sample_Report.csv";
}

my ($projectIndex,$fileIndex,$typeIndex,$sampleNameIndex,$genomeIndex,$laneIndex);
my ($project,$type,$sample,$genome,$fastq,$lane);
my (%samples,%genomes,%types,%files,%projects);
open(IN,$sampleSheet);
while(<IN>){
    chomp $_;
    if ($_ =~ /^output/){
	my @labels = split(/\,/,$_);
	for (my $i=0;$i<=$#labels;$i++){
	    if ($labels[$i] eq "order"){
		$projectIndex = $i;
	        print STDERR "Project name in Column $i\n";
	    }
	    elsif ($labels[$i] eq "output"){
		$fileIndex = $i;
		print STDERR "fastq in Column $i\n";
	    }
	    elsif ($labels[$i] eq "lane"){
		$laneIndex = $i;
		print STDERR "lane in Column $i\n";
	    }
	    elsif ($labels[$i] eq "order type"){
		$typeIndex = $i;
		print STDERR "Order type in Column $i\n";
	    }
	    elsif (($labels[$i] eq "sample name") || ($labels[$i] eq "sample.name")){
		$sampleNameIndex = $i;
		print STDERR "Sample name in Column $i\n";
	    }
	    elsif ($labels[$i] eq "reference"){
		print STDERR "Reference name in Column $i\n";
		$genomeIndex = $i;
	    }
	}
    }else {
	$type = "default";
	my @data = split(/\,/,$_);
	if ($projectIndex){
	    $project = $data[$projectIndex];
	} else { $project = $projectID;}
	if ($typeIndex){
	    $type = $data[$typeIndex];
	    if ($type eq "ChIP-Seq"){ $type = "chipseq";}
	    elsif ($type =~ /RNA/){ $type = "RNA";}
	}
	$sample = $data[$sampleNameIndex];
	$lane = $data[$laneIndex];
	$fastq = $data[$fileIndex];
	if ($genomeIndex){
	    $genome = $data[$genomeIndex];
	}
	$genomes{$project} = $genome;
	if ($genome ne ""){
	    $sample = $sample."\_S$lane\_$genome\_$type";
	}else {
	    $sample = $sample."\_S$lane";
	}
	if (!exists($samples{$project})){
	    $samples{$project} = $sample;
	} else {
	    $samples{$project} .= ",$sample";
	}
	$files{$sample} = $fastq;
	$types{$project} = $type;
#	print "$project\t$sample\t$fastq\t$type\t$genome\n";
    }
}

foreach my $project (keys(%samples)){
    if ($projectID ne ""){
	if ($project eq $projectID){
	    print STDERR "$project\t$types{$project}\t$genomes{$project}\t$samples{$project}\n\n";
	    `mkdir $outputDir\/$project`;
	    my @samples = split(/\,/,$samples{$project});
	    foreach my $sample (@samples){
		my $cmd = "ln -s $fastQdir\/$files{$sample} $outputDir\/$project\/$sample.fastq.gz";
#		print STDERR "$cmd\n";
		system($cmd);
	    }
	}
    } else{
	# Just do all Project IDs
	print STDERR "$project\t$types{$project}\t$genomes{$project}\t$samples{$project}\n\n";
	`mkdir $outputDir\/$project`;
	my @samples = split(/\,/,$samples{$project});
	foreach my $sample (@samples){
	    if (-e "$fastQdir\/$files{$sample}"){
		my $cmd = "ln -s $fastQdir\/$files{$sample} $outputDir\/$project\/$sample.fastq.gz";
#		print STDERR "$cmd\n";
		system($cmd);
	    } else {
		print STDERR "ERR:  Could not find $fastQdir\/$files{$sample}! Skipping it.\n";
	    }
	}
    }
}
