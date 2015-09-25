#!/bin/bash
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
#
# This script is for identify the regions where there are significant changes between the two libraries.
# The two libraries can be KO vs WT, or before or after differentiation or  normal and disease states.1
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  11/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

SICER=/projects/b1025/tools/SICER_V1.1/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 8 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["KO bed file"] ["KO control file"] ["WT bed file"]  ["WT control file"] ["window size (bp)"] ["gap size (bp)"] ["FDR for KO vs KOCONTROL or WT vs WTCONTROL"] ["FDR for WT vs KO"] 
    echo ""
    exit 1
fi

TEMP=`expr $6 % $5`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
echo "#############################################"
echo "######           SICER v1.1            ######"
echo "#############################################"
# Input KO bed file
KOBED=$1
KO=${KOBED%.*}
echo "KO library: $KOBED"


KOCONTROLBED=$2
KOCONTROL=${KOCONTROLBED%.*}
echo "KO control library: $KOCONTROLBED"


# Input WT bed file.
WTBED=$3
WT=${WTBED%.*}
echo "WT library: $WTBED"


# Input WT bed file.
WTCONTROLBED=$4
WTCONTROL=${WTCONTROLBED%.*}
echo "WT library: $WTCONTROLBED"



# Species, for allowed species see GenomeData.py
SPECIES=hg18
echo "Species: $SPECIES"

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=0.74
echo "Effective genome size as a fraction of reference genome: $EFFECTIVEGENOME"

# KOTHRESHOLD is the threshold is for redundancy allowed for reads. KOTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
KOTHRESHOLD=1
echo "Threshold for redundancy allowed for treated reads: $KOTHRESHOLD"

# WTTHRESHOLD is the threshold is for redundancy allowed for reads. WTTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
WTTHRESHOLD=1
echo "Threshold for redundancy allowed for WT reads: $WTTHRESHOLD"

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$5
echo "Window size: $WINDOW_SIZE bps"

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=150
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"

#GAP_SIZE is in base pairs.
GAP_SIZE=$6
echo "Gap size: $GAP_SIZE bps"

#EVALUE is the number of islands expected in random background. The E value is used for identification of candidate islands that exhibit clustering.
FDR=$7
echo "FDR for identification of enriched islands: $EVALUE"

#False discovery rate WTling significance
FDR_WT_KO=$8
echo "FDR for identification of significant changes: $FDR_WT_KO "
#================================================================================
#############################################
# ######  SET UP DIRECTORY STRUCTURE ###### #
#############################################
# If data files are not in the current directory, replace this with
# the path to the data files.
DATADIR=.
KODIR=$DATADIR
WTDIR=$DATADIR
# If You want the output files not in the current directory, replace this.
OUTPUTDIR=.
KOOUTPUTDIR=$OUTPUTDIR
WTOUTPUTDIR=$OUTPUTDIR


#================================================================================
#############################################
# ###### DEFINITION OF OUTPUT FILES  ###### #
#############################################


# This file stores the preprocessed raw bed file.
FILTEREDKOBED=$KO-$KOTHRESHOLD-removed.bed
FILTEREDWTBED=$WT-$WTTHRESHOLD-removed.bed

#This file stores the candidate islands in  "chr     start   end     ChIP_island_read_count" format
KOISLAND=$KO-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-island.bed
WTISLAND=$WT-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-island.bed

UNIONISLAND=$KO-vs-$WT-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-union.island

# This file stores the island-filtered non-redundant raw reads 
KOISLANDFILTEREDRAWBED=$KO-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.bed
WTISLANDFILTEREDRAWBED=$WT-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.bed

#This file stores the summary of candidate islands, including tags counts, 
# pvalue, fold change and BH-corrected p-value
MERGEDISLANDSUMMARYFILE=$KO-and-$WT-W$WINDOW_SIZE-G$GAP_SIZE-summary

#This file stores the summary of significant islands identified with BH-corrected p-value criterion.
INCREASEDISLANDS=$KO-W$WINDOW_SIZE-G$GAP_SIZE-increased-islands-summary-FDR$FDR_WT_KO
DECREASEDISLANDS=$KO-W$WINDOW_SIZE-G$GAP_SIZE-decreased-islands-summary-FDR$FDR_WT_KO


echo " "
echo " "
echo " "
echo " "
echo "Process the $KO and $KOCONTROLBED libraries using SICER.sh"
echo "sh $SICER/SICER.sh $KODIR $KOBED $KOCONTROLBED $KOOUTPUTDIR  $SPECIES $KOTHRESHOLD $WINDOW_SIZE $FRAGMENT_SIZE $EFFECTIVEGENOME $GAP_SIZE $FDR"
sh $SICER/SICER.sh $KODIR $KOBED $KOCONTROLBED $KOOUTPUTDIR  $SPECIES $KOTHRESHOLD $WINDOW_SIZE $FRAGMENT_SIZE $EFFECTIVEGENOME $GAP_SIZE $FDR


echo " "
echo " "
echo " "
echo " "
echo "Process the $WT and $WTCONTROLBED libraries using SICER.sh"
echo "sh $SICER/SICER.sh $WTDIR $WTBED $WTCONTROLBED $WTOUTPUTDIR $SPECIES $WTTHRESHOLD $WINDOW_SIZE $FRAGMENT_SIZE $EFFECTIVEGENOME $GAP_SIZE $FDR"
sh $SICER/SICER.sh $WTDIR $WTBED $WTCONTROLBED $WTOUTPUTDIR $SPECIES $WTTHRESHOLD $WINDOW_SIZE $FRAGMENT_SIZE $EFFECTIVEGENOME $GAP_SIZE $FDR


echo " "
echo " "
echo ""
echo ""
echo "Merge the two identified sets of significant islands..."
echo "python $SICER/src/find_union_islands.py -s $SPECIES -a $WTOUTPUTDIR/$WTISLAND -b $KOOUTPUTDIR/$KOISLAND -o $OUTPUTDIR/$UNIONISLAND"
python $SICER/src/find_union_islands.py -s $SPECIES -a $WTOUTPUTDIR/$WTISLAND -b $KOOUTPUTDIR/$KOISLAND -o $OUTPUTDIR/$UNIONISLAND


echo ""
echo ""
echo "Find island-filtered read counts from the two libraries on the merged islands  and calculate significance of increase and decrease ..."

#Use preprocessed raw reads for comparision
echo "python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES  -a $KODIR/$FILTEREDKOBED -b $WTDIR/$FILTEREDWTBED -d $OUTPUTDIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE"
python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES -a $KODIR/$FILTEREDKOBED -b $WTDIR/$FILTEREDWTBED -d $OUTPUTDIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE

#Use island filtered reads for comparison
#echo "python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES  -a $KOOUTPUTDIR/$KOISLANDFILTEREDRAWBED -b $WTOUTPUTDIR/$WTISLANDFILTEREDRAWBED -d $OUTPUTDIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE"
#python $SICER/src/compare_two_libraries_on_islands.py -s $SPECIES  -a $KOOUTPUTDIR/$KOISLANDFILTEREDRAWBED -b $WTOUTPUTDIR/$WTISLANDFILTEREDRAWBED -d $OUTPUTDIR/$UNIONISLAND -f $FRAGMENT_SIZE -o $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE


echo ""
echo ""
echo "Identify significantly increased islands using BH corrected p-value cutoff ..."
echo "python $SICER/src/filter_islands_by_significance.py -i $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE -p $FDR_WT_KO -c 9 -o  $OUTPUTDIR/$INCREASEDISLANDS"
python $SICER/src/filter_islands_by_significance.py -i $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE -p $FDR_WT_KO -c 9 -o  $OUTPUTDIR/$INCREASEDISLANDS


echo ""
echo ""
echo "Identify significantly decreased islands using BH-corrected p-value cutoff ..."
echo "python $SICER/src/filter_islands_by_significance.py -i $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE -p $FDR_WT_KO -c 12 -o  $OUTPUTDIR/$DECREASEDISLANDS"
python $SICER/src/filter_islands_by_significance.py -i $OUTPUTDIR/$MERGEDISLANDSUMMARYFILE -p $FDR_WT_KO -c 12 -o  $OUTPUTDIR/$DECREASEDISLANDS



echo "Done!"
