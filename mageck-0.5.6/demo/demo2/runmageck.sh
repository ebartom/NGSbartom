#!/bin/bash


# The old "run" command will be discontinued.
# mageck run --fastq test1.fastq test2.fastq -l library.txt -n demo --sample-label L1,CTRL -t L1 -c CTRL

# instead, execute "count" command first to generate read count tables
mageck count -l library.txt -n demo --sample-label L1,CTRL  --fastq test1.fastq test2.fastq 

# then, run "test" command 
mageck test -k demo.count.txt -t L1 -c CTRL -n demo

