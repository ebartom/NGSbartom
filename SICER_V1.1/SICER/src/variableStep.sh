#!/bin/sh


#SPAN=200
#$1 Input file name
S=$1
#$2 Output file name
R=$2
#$3 Track label
#$4 window size
SPAN=$4

rm -f ${R}
head -45 ${S} | egrep "^browser|^track" > ${R}
echo "track type=wiggle_0 name=$3" >> ${R}
grep "^chr" ${S} | cut -f1 | sort -u > chr.list
cat chr.list | while read C
do
    echo "variableStep chrom=${C} span=${SPAN}" >> ${R}
    awk '{if (match($1,"^'"${C}"'$")) { print } }' ${S} | sort -k2n | awk '
{
    printf "%d\t%g\n", $2+1, $4
}
' >> ${R}
done
