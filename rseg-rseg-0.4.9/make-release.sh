#! env bash 

# make-release.sh 
# Copyright 2012 Song Qiang <qiang.song@usc.edu>
#
# pull auxillary library files and setting up Makefiles

RSEG_ROOT=$PWD

cd $RSEG_ROOT/src
svn export https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools
svn export http://smithlab.usc.edu/repos/smithlab_cpp

cd $RSEG_ROOT/src/rseg
cat Makefile |sed 's|^ifndef SMITHLAB_CPP$|SMITHLAB_CPP = ../smithlab_cpp\n&|' \
	> tmp && mv tmp Makefile

cd $RSEG_ROOT