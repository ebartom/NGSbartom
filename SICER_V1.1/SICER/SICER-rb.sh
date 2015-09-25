#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  10/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################
PATHTO=/projects/b1025/tools/SICER_V1.1
SICER=$PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 10 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["InputDir"] ["bed file"] ["OutputDir"] ["species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["E-value"] 
    echo ""
    exit 1
fi

#TEMP=$[$3 % $2] 
TEMP=`expr $9 % $6`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
#Parameters for running SICER without control library and using random background model.

# Input directory
DATADIR=$1
DATADIR=${DATADIR:=.}
SAMPLEDIR=$DATADIR

# Input sample bed file
SAMPLEBED=$2
SAMPLE=${SAMPLEBED%.*}

#Output directory
OUTPUTDIR=$3
OUTPUTDIR=${OUTPUTDIR:=.}
SUMMARY_DIR=$OUTPUTDIR
ISLANDS_BED_DIR=$OUTPUTDIR

# Species, for allowed species see GenomeData.py
SPECIES=$4
SPECIES=${SPECIES:=hg18}

# THRESHOLD is the threshold is for redundancy for reads. THRESHOLD=n
# means that each read has at most n copy after preprocessing.
THRESHOLD=$5
THRESHOLD=${THRESHOLD:=1}

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$6

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=$7
FRAGMENT_SIZE=${FRAGMENT_SIZE:=150}

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=$8
EFFECTIVEGENOME=${EFFECTIVEGENOME:=0.74}

if [ $(echo "$EFFECTIVEGENOME > 1"|bc) -eq 1 ]; then
	echo ""
	echo " $EFFECTIVEGENOME needs to be between 0 and 1 "
	echo "" 
	exit 1
fi

#GAP_SIZE is in base pairs.
GAP_SIZE=$9

#EVALUE is the number of islands expected in random background. The
#EVALUE is used to determine the score threshold s_T.
EVALUE=${10}


echo "#############################################"
echo "######           SICER v1.1            ######"
echo "#############################################"

echo "Input library directory: $DATADIR"
echo "ChIP library: $SAMPLEBED"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Threshold for redundancy allowed for reads: $THRESHOLD"
echo "Window size: $WINDOW_SIZE bps"
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"
echo "Effective genome size as a fraction of the reference genome of $SPECIES: $EFFECTIVEGENOME"
echo "Gap size: $GAP_SIZE bps"
echo "Evalue for identification of significant islands: $EVALUE"
#================================================================================


# This file stores the preprocessed raw bed file.
FILTEREDSAMPLEBED=$SAMPLE-${THRESHOLD}-removed.bed

# This file stores the summary graph.  
SUMMARY=$SAMPLE-W$WINDOW_SIZE.graph
NORMALIZEDSUMMARY=$SAMPLE-W$WINDOW_SIZE-normalized.graph
# This file stores the histogram of window read-count.  
SUMMARYHIST=$SAMPLE-W$WINDOW_SIZE.graphhist
# This file stores the summary graph in wig vstep format
NORMALIZEDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-normalized.wig

#This file stores the candidate islands.
ISLAND=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE.scoreisland
# This file stores the windows on islands. 
FILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-islandfiltered.scoregraph
# This file stores the histogram of island scores
ISLANDSCOREHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE.islandscorehist
# This file stores the histogram of island lengths
ISLANDLENGTHHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE.islandlengthhist

# This file stores the island-filtered non-redundant raw reads 
ISLANDFILTEREDRAWBED=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-islandfiltered.bed
# This file stores summary graph made by the island-filtered non-redundant raw reads 
ISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-islandfiltered.graph
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads 
NORMALIZEDISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-islandfiltered-normalized.graph 
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads in wig vstep format
NORMALIZEDISLANDFILTEREDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-E$EVALUE-islandfiltered-normalized.wig



echo " "
echo " "
echo "Preprocess the raw $SAMPLE file to remove redundancy with threshold $THRESHOLD..."
echo "python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $THRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED"
python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $THRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED

echo " "
echo " "
echo "Partion the genome in windows ..."
echo "Generate summary files for $SAMPLE..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SUMMARY 

echo ""
echo ""
echo "Normalize summary graph by total redundancy-removed reads per million for $SAMPLE ..."
echo "python $SICER/src/normalize.py -i $SUMMARY_DIR/$SUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSUMMARY"
python $SICER/src/normalize.py -i $SUMMARY_DIR/$SUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSUMMARY

echo ""
echo ""
echo "Convert the normalized summary graph into wig vstep format..."
echo "sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSUMMARY $SUMMARY_DIR/$NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSUMMARY $SUMMARY_DIR/$NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE

rm $SUMMARY_DIR/$NORMALIZEDSUMMARY

## Generating window read-count statistics, optional
#     echo " "
#     echo " "
#     echo "Generate window read count histograms for $SAMPLE ... "
#     echo "python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$SUMMARYHIST"
#     python $SICER/utility/get_windows_histogram.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -o $SUMMARY_DIR/$SUMMARYHIST



echo " "
echo " "
echo "Find significant islands with E-value $EVALUE for $SAMPLE..."
echo "python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE -f $ISLANDS_BED_DIR/$ISLAND"
python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE  -f $ISLANDS_BED_DIR/$ISLAND

# Generating island statistics, optional
# echo " "
# echo " "
# echo "Get island statistics for $SAMPLE..."
# echo "python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -i $ISLANDS_BED_DIR/$ISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$ISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$ISLANDLENGTHHIST -r None"
# python $SICER/utility/islands_statistics_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -i $ISLANDS_BED_DIR/$ISLAND -w $WINDOW_SIZE -o $ISLANDS_BED_DIR/$ISLANDSCOREHIST   -q $ISLANDS_BED_DIR/$ISLANDLENGTHHIST -r None

echo ""
echo ""
echo "Filter reads with identified significant islands for $SAMPLE..."
echo "python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLAND  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED"
python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED  -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLAND  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED

echo ""
echo ""
echo "Make summary graph with filtered reads for $SAMPLE..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY


echo ""
echo ""
echo "Normalize summary graph with filtered reads for $SAMPLE by total island filtered reads per million..."
echo "python $SICER/src/normalize.py -i $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY -a 3 -t 1000000 -o $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY"
python $SICER/src/normalize.py -i $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY -a 3 -t 1000000 -o $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY

echo ""
echo ""
echo "Convert the summary graph made with the filtered reads into wig vstep format and normalize by total island-filtered read count per million..."
echo "sh $SICER/src/variableStep.sh $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARYWIG $SAMPLE-islandfiltered $WINDOW_SIZE

rm $ISLANDS_BED_DIR/$ISLANDFILTEREDSUMMARY
rm $ISLANDS_BED_DIR/$NORMALIZEDISLANDFILTEREDSUMMARY

echo ""
echo ""
echo "Done!"
