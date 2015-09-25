#!/bin/bash

# !!!!Replace it with your own directory!!!
SICER=/home/data/SICER1.1/SICER

# ["InputDir"] ["bed file"] ["OutputDir"] ["species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["E-value"]
sh $SICER/SICER-rb.sh /home/data/SICER1.1/SICER/ex test.bed . hg18 1 200 150 0.74 400 100