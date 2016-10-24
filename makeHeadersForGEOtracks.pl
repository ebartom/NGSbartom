#!/usr/bin/perl -w
use strict;
use File::Basename;

my $trackDir = $ARGV[0];
my $urlBase = $ARGV[1];

if ($#ARGV <1){
    die "Usage:  $0 <trackDirectory> <baseURL>\n";
}

# As an example:
# perl makeHeadersForGEOtracks.pl /projects/b1025/etb/deqingMLL2/onGEO/ https://s3-us-west-2.amazonaws.com/ash-tracks/etb/etb.MLL2GEO/

my @trackfiles = split(/\n/,`ls $trackDir/*.bw`);

foreach my $bw (@trackfiles){
    my $name = fileparse($bw);
    $name =~ s/.bw$//g;
    $urlBase =~ s/\/$//g;
    print "track type=bigWig name=$name.bw description=$name.rpm graphtype=bar maxHeightPixels=128:60:11 visibility=full color=0,0,255 itemRGB=on autoScale=on bigDataUrl=$urlBase/$name.bw\n";
}
