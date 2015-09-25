#!/bin/bash
#MSUB -l nodes=1:ppn=1
#MSUB -A b1025
#MSUB –m abe
#MSUB -j oe
#MSUB –l walltime=24:00:00
#MSUB -q normal

export PATH=$PATH:/projects/b1025/tools/bedtools2/bin/

Files="Trr-C2398A-H3K4me1-150505
Trr-C2398A-Input-150505
Trr-WT-H3K4me1-150505
Trr-WT-Input-150505
Trr-Y2383F-H3K4me1-150505
Trr-Y2383F-Input-150505
"
bamDir="/projects/b1025/tango/TANGO-067.077.078/TANGO-077/bam"
outDir="/projects/b1025/arw/test/chip/sicer/bed"

[ ! -d ${outDir} ] && mkdir -p ${outDir}

for File in $Files
do
    bedtools bamtobed -i ${bamDir}/${File}.bam >> ${outDir}/${File}.bed
done
