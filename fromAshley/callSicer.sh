#!/bin/bash
#MSUB -l nodes=1:ppn=1
#MSUB -A b1025
#MSUB –m abe
#MSUB -j oe
#MSUB –l walltime=24:00:00
#MSUB -q normal

module load python
export PATH=$PATH:/projects/b1025/tools/SICER_V1.1/SICER/

ref="dm3"
orderLog="_077_log"
bedDir="/projects/b1025/arw/test/chip/sicer/bed"
outDir="/projects/b1025/arw/test/chip/sicer"

Files="Trr-C2398A
Trr-WT
Trr-Y2383F
"

for File in $Files
do
SICER.sh ${bedDir} ${File}-H3K4me1-150505.bed ${File}-Input-150505.bed ${outDir} ${ref} 1 200 150 0.80 600 1e-8 &> ${outDir}/${File}${orderLog}
done
