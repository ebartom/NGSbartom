#!/bin/bash

# This script is a wrapper to run buildPipelineScripts.pl, accepting and passing only the "build" parameters.

# The following parameters are supported:

# -o <outputDirectory>
# -bs <baseSpaceDirectory>
# -f <fastqDirectory>
# -bam <bamDirectory>
# -c <comparisons.csv file>
# -4C <4C description file>
# -chip <ChIP Description file>
# -p <numProcessors>
# -m <1|0>
# -a <aligner>
# -g <assembly/genome>
# -t <RNA|chipseq|4C>
# -buildBcl2fq <1|0>
# -runTrim <1|0>
# -buildAlign <1|0>
# -makeTracks <1|0>
# -buildEdgeR <1|0>
# -buildPeakCaller <1|0>
# -buildDiffPeaks <1|0>
# -build4C <1|0>

options=$(/usr/bin/getopt -a -o o:f:c:p:m:a:g:t: -l bs:,bam:,4C:,chip:,buildBcl2fq:,runTrim:,buildAlign:,makeTracks:,buildEdgeR:,buildPeakCaller:,buildDiffPeaks:,build4C: -- "$@")

if [ $? -ne 0 ] 
then
    echo "Incorrect options provided."
    exit 1
fi

# need to find output and basespace directories
eval set -- "$options"
while true; do
	case "$1" in
	-o)
		OUTPUT_DIR=$2; shift; shift;;
	-bs)
		BASESPACE_DIR=$2; shift; shift;;
	--)
		shift; break;;
	*)
		shift;;
	esac
done

echo "OUTPUT_DIR set to $OUTPUT_DIR"

echo "BASESPACE_DIR set to $BASESPACE_DIR"

# have to remove quotes from arguments
options=$( echo $options | sed "s/'//g" )

# make sure module system is loaded
if [ "$(type -t module)" != "function" ]
then
	MODULEPATH=""
	. /etc/profile
fi

echo "calling buildPipelineScripts.pl with arguments: -uploadASHtracks 0 $options"

./buildPipelineScripts.pl -uploadASHtracks 0 $options

echo "Making all scripts executable."
find ${OUTPUT_DIR}/*/scripts -iname '*.sh' -exec chmod 755 {} +

date
echo "Running ${baseSpaceDirectory}/runBcl2fq.sh if exists"
if [ -f $BASESPACE_DIR/runBcl2fq.sh ]
then
	/bin/bash -l $BASESPACE_DIR/runBcl2fq.sh
fi

date
echo "run $OUTPUT_DIR/$sample_project/scripts/run_4C_demultiplex.sh if exists"
find ${OUTPUT_DIR}/*/scripts -name run_4C_demultiplex.sh -exec /bin/bash -l {} \;

date
echo "run alignment scripts"
find ${OUTPUT_DIR}/*/scripts -iname '*align.sh' -print | xargs -L 1 -P $(find ${OUTPUT_DIR}/*/scripts -iname '*align.sh' -print | wc -l) /bin/bash -l

date
echo "run genotype scripts"
find ${OUTPUT_DIR}/*/scripts -iname '*genotype.sh' -print | xargs -L 1 -P $(find ${OUTPUT_DIR}/*/scripts -iname '*genotype.sh' -print | wc -l) /bin/bash -l

date
echo "run downstream rna analysis"
find ${OUTPUT_DIR}/*/scripts -iname 'downstreamRNAanalysis.sh' -exec /bin/bash -l {} \;

date
echo "run peak caller scripts"
find ${OUTPUT_DIR}/*/scripts -name '*callPeaks.sh' -print | xargs -L 1 -P $(find ${OUTPUT_DIR}/*/scripts -name '*callPeaks.sh' -print | wc -l) /bin/bash -l

date
echo "run diffPeak scripts"
find ${OUTPUT_DIR}/*/scripts -name '*diffPeaks.sh ' -print | xargs -L 1 -P $(find ${OUTPUT_DIR}/*/scripts -name '*diffPeaks.sh' -print | wc -l) /bin/bash -l

date
echo "Done"
