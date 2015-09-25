#!/bin/bash

# !!!!Replace it with your own directory!!!
SICER=/home/data/SICER1.1/SICER

#["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]
sh $SICER/SICER.sh /home/data/SICER1.1/SICER/ex test.bed control.bed . hg18 1 200 150 0.74 600 .01