#!/bin/bash
# 
# Authors: Chongzhi Zang, Weiqun Peng
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
#
# Version 1.1  6/9/2010

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################

PATHTO=/projects/b1025/tools/SICER_V1.1
SICER=$PATHTO/SICER
PYTHONPATH=$SICER/lib
export PYTHONPATH

if [ $# -lt 11 ]; then
    echo ""
    echo 1>&2 Usage: $0 ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]  
    echo ""
    exit 1
fi


TEMP=`expr ${10} % $7`

if [ $TEMP != 0 ]; then
	echo ""
	echo "Gap_size needs to be multiples of window_size, ie, 0, 2*window_size, etc "
	echo ""
	exit 1
fi

#================================================================================
#############################################
# ######  SET UP DIRECTORY STRUCTURE ###### #
#############################################

# The path to the data files.
DATADIR=$1
SAMPLEDIR=$DATADIR
CONTROLDIR=$DATADIR

# Output directory
OUTPUTDIR=$4
SUMMARY_DIR=$OUTPUTDIR
ISLANDS_BED_DIR=$OUTPUTDIR
#================================================================================

# Input sample bed file
SAMPLEBED=$2
SAMPLE=${SAMPLEBED%.*}


# Input control bed file
CONTROLBED=$3
CONTROL=${CONTROLBED%.*}



# Species, for allowed species see GenomeData.py
SPECIES=$5
SPECIES=${SPECIES:=hg18}

# CHIPTHRESHOLD is the threshold is for redundancy allowed for reads. CHIPTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
CHIPTHRESHOLD=$6
CHIPTHRESHOLD=${CHIPTHRESHOLD:=1}


# CONTROLTHRESHOLD is the threshold is for redundancy allowed for reads. CONTROLTHRESHOLD=n
# means that each read has at most n copy after preprocessing.
CONTROLTHRESHOLD=$6
CONTROLTHRESHOLD=${CONTROLTHRESHOLD:=1}

# WINDOW_SIZE is the size of the windows to scan the genome width. 
# One WINDOW_SIZE is the smallest possible island.
WINDOW_SIZE=$7

# FRAGMENT_SIZE is for determination of the amound of shift for center
# of the DNA fragments represented by the reads. FRAGMENT_SIZE=150
# means the shift is 75.
FRAGMENT_SIZE=$8
FRAGMENT_SIZE=${FRAGMENT_SIZE:=150}

# Effective Genome as fraction of the genome size. It depends on read length.
EFFECTIVEGENOME=$9
EFFECTIVEGENOME=${EFFECTIVEGENOME:=0.74}

#EVALUE is the number of islands expected in random background. The E value is used for identification of candidate islands that exhibit clustering.
EVALUE=1000

#GAP_SIZE is in base pairs.
GAP_SIZE=${10}

#False discovery rate controlling significance
FDR=${11}


echo "#############################################"
echo "######           SICER v1.1            ######"
echo "#############################################"

echo "Input library directory: $DATADIR"
echo "ChIP library: $SAMPLEBED"
echo "Control library: $CONTROLBED"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Threshold for redundancy allowed for chip reads: $CHIPTHRESHOLD"
echo "Threshold for redundancy allowed for control reads: $CONTROLTHRESHOLD"
echo "Window size: $WINDOW_SIZE bps"
echo "Fragment size: $FRAGMENT_SIZE bps. The shift for reads is half of $FRAGMENT_SIZE"
echo "Effective genome size as a fraction of the reference genome of $SPECIES: $EFFECTIVEGENOME"
echo "Gap size: $GAP_SIZE bps"
echo "Evalue for identification of candidate islands that exhibit clustering: $EVALUE"
echo "False discovery rate controlling significance: $FDR"

#================================================================================

#############################################
# ###### DEFINITION OF OUTPUT FILES  ###### #
#############################################

# This file stores the preprocessed raw bed file.
FILTEREDSAMPLEBED=$SAMPLE-${CHIPTHRESHOLD}-removed.bed
FILTEREDCONTROLBED=$CONTROL-${CONTROLTHRESHOLD}-removed.bed

# This file stores the summary graph.  
SUMMARY=$SAMPLE-W$WINDOW_SIZE.graph
NORMALIZEDSUMMARY=$SAMPLE-W$WINDOW_SIZE-normalized.graph
# This file stores the histogram of window read-count.  
SUMMARYHIST=$SAMPLE-W$WINDOW_SIZE.graphhist
# This file stores the summary graph in wig vstep format
NORMALIZEDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-normalized.wig

#This file stores the candidate islands.
ISLAND=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE.scoreisland

#This file stores the summary of candidate islands, including chrom start end read-count_sample read-count-control pvalue, fold change and qvalue
ISLANDSIGNIFICANCEFILE=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-islands-summary

#This file stores the summary of significant islands identified with FDR criterion.
SIGNIFICANTISLANDS=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-islands-summary-FDR$FDR


#This file stores the significant islands in  "chr     start   end" format
ISLANDBED=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-island.bed
# This file stores the histogram of island scores
ISLANDSCOREHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR.islandscorehist
# This file stores the histogram of island lengths
ISLANDLENGTHHIST=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR.islandlengthhist


# This file stores the island-filtered non-redundant raw reads 
ISLANDFILTEREDRAWBED=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.bed
# This file stores summary graph made by the island-filtered non-redundant raw reads 
ISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered.graph
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads 
NORMALIZEDISLANDFILTEREDSUMMARY=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered-normalized.graph 
# This file stores normalized summary graph made by the island-filtered non-redundant raw reads in wig vstep format
NORMALIZEDISLANDFILTEREDSUMMARYWIG=$SAMPLE-W$WINDOW_SIZE-G$GAP_SIZE-FDR$FDR-islandfiltered-normalized.wig

echo " "
echo " "
echo "Preprocess the raw $SAMPLE file to remove redundancy with threshold $CHIPTHRESHOLD..."
echo "python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $CHIPTHRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED"
python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $SAMPLEDIR/$SAMPLEBED -t $CHIPTHRESHOLD -o $OUTPUTDIR/$FILTEREDSAMPLEBED

echo " "
echo " "
echo "Preprocess the raw $CONTROL file to remove redundancy with threshold $CONTROLTHRESHOLD..."
echo "python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $CONTROLDIR/$CONTROL.bed -t $CONTROLTHRESHOLD -o $OUTPUTDIR/$FILTEREDCONTROLBED"
python $SICER/src/remove_redundant_reads.py -s $SPECIES -b $CONTROLDIR/$CONTROL.bed -t $CONTROLTHRESHOLD -o $OUTPUTDIR/$FILTEREDCONTROLBED


echo " "
echo " "
echo "Partion the genome in windows ..."
echo "Generate summary files ..."
echo "python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SUMMARY"
python $SICER/src/run-make-graph-file-by-chrom.py -s $SPECIES -b $OUTPUTDIR/$FILTEREDSAMPLEBED -w $WINDOW_SIZE -i $FRAGMENT_SIZE -o $SUMMARY_DIR/$SUMMARY 

echo ""
echo ""
echo "Normalize summary graph by total island filtered reads per million for $SAMPLE ..."
echo "python $SICER/src/normalize.py -i $SUMMARY_DIR/$SUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSUMMARY"
python $SICER/src/normalize.py -i $SUMMARY_DIR/$SUMMARY -a 3 -t 1000000 -o $SUMMARY_DIR/$NORMALIZEDSUMMARY

echo ""
echo ""
echo "Convert the normalized summary graph into wig vstep format..."
echo "sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSUMMARY $SUMMARY_DIR/$NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE"
sh $SICER/src/variableStep.sh $SUMMARY_DIR/$NORMALIZEDSUMMARY $SUMMARY_DIR/$NORMALIZEDSUMMARYWIG $SAMPLE $WINDOW_SIZE

rm $SUMMARY_DIR/$NORMALIZEDSUMMARY


echo " "
echo " "
echo "Find candidate islands exhibiting clustering ..."
echo "python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE -f $ISLANDS_BED_DIR/$ISLAND"
python $SICER/src/find_islands_in_pr.py -s $SPECIES -b $SUMMARY_DIR/$SUMMARY -w $WINDOW_SIZE -g $GAP_SIZE -t $EFFECTIVEGENOME -e $EVALUE  -f $ISLANDS_BED_DIR/$ISLAND



echo ""
echo ""
echo "Calculate significance of candidate islands using the control library ..."
echo "python $SICER/src/associate_tags_with_chip_and_control_w_fc_q.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$ISLAND -f $FRAGMENT_SIZE -t $EFFECTIVEGENOME -o $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE"
python $SICER/src/associate_tags_with_chip_and_control_w_fc_q.py -s $SPECIES  -a $OUTPUTDIR/$FILTEREDSAMPLEBED -b $OUTPUTDIR/$FILTEREDCONTROLBED -d $ISLANDS_BED_DIR/$ISLAND -f $FRAGMENT_SIZE -t $EFFECTIVEGENOME -o $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE

echo ""
echo ""
echo "Identify significant islands using FDR criterion ..."
#echo "python $SICER/src/find_significant_islands.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -q $FDR -o  $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS"
#python $SICER/src/find_significant_islands.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -q $FDR -o  $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS
echo "python $SICER/src/filter_islands_by_significance.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -p $FDR -c 7 -o $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS"
python $SICER/src/filter_islands_by_significance.py -i $ISLANDS_BED_DIR/$ISLANDSIGNIFICANCEFILE -p $FDR -c 7 -o $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS



echo ""
echo ""
echo "Convert island summary to island bed file of format chr start end ChIP-read-count"
echo "python $SICER/utility/convert_summary_to_bed.py -i $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS  -o  $ISLANDS_BED_DIR/$ISLANDBED"
python $SICER/utility/convert_summary_to_bed.py -i $ISLANDS_BED_DIR/$SIGNIFICANTISLANDS  -o  $ISLANDS_BED_DIR/$ISLANDBED


echo ""
echo ""
echo "Filter reads with identified significant islands..."
echo "python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLANDBED  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED"
python $SICER/utility/filter_raw_tags_by_islands.py -s $SPECIES -a $OUTPUTDIR/$FILTEREDSAMPLEBED -i $FRAGMENT_SIZE -b $ISLANDS_BED_DIR/$ISLANDBED  -o  $ISLANDS_BED_DIR/$ISLANDFILTEREDRAWBED


echo ""
echo ""
echo "Make summary graph with filtered reads..."
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
