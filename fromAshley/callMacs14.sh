#!/bin/bash
#MSUB -l nodes=1:ppn=1
#MSUB -A b1025
#MSUB –m abe
#MSUB -j oe
#MSUB –l walltime=24:00:00
#MSUB -q normal

export PATH=$PATH:/projects/b1025/tools/MACS-1.4.2/bin
export PYTHONPATH=/projects/b1025/tools/MACS-1.4.2/lib/python2.6/site-packages:$PYTHONPATH

Files="Trr-C2398A
Trr-WT
Trr-Y2383F
"

ref="dm"
orderLog="_077_log"
bamDir="/projects/b1025/tango/TANGO-067.077.078/TANGO-077/bam"
outDir="/projects/b1025/arw/test/chip/macs14"

for File in $Files
do
    macs14 -t ${bamDir}/${File}-H3K4me3-150505.bam -c ${bamDir}/${File}-Input-150505.bam -f BAM -g ${ref} -n ${outDir}/${File}H3K4me3_150505_077 2> ${outDir}/${File}H3K4me3_150505${orderLog}
done

