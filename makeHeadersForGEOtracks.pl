#!/usr/bin/perl -w
use strict;
use File::Basename;

my $trackDir = $ARGV[0];
my $urlBase = $ARGV[1];
my $uploadTracks = $ARGV[2];

if ($#ARGV <1){
    die "Usage:  $0 <trackDirectory> <baseURL> <uploadTracks>\n";
}

# As an example:
# perl makeHeadersForGEOtracks.pl /projects/b1025/etb/deqingMLL2/onGEO/ https://s3-us-west-2.amazonaws.com/ash-tracks/etb/etb.MLL2GEO/ 1

my @trackfiles = split(/\n/,`ls $trackDir/*.bw`);
my $suffix = "ash-tracks/TANGO/";
if ($urlBase =~ /amazonaws.com\/([\w\-\/\.\_]+)$/){
    $suffix = $1;
}

foreach my $bw (@trackfiles){
    my $name = fileparse($bw);
    $name =~ s/.bw$//g;
    $urlBase =~ s/\/$//g;
    if ($uploadTracks == 1){
	my $returned = system("aws s3 cp $bw s3://$suffix --region us-west-2");
    }
    print "track type=bigWig name=$name.bw description=$name.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=$urlBase/$name.bw\n";
}
