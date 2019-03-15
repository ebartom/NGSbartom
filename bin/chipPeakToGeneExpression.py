#!/software/anaconda3.6/bin/python3
# Author: Patrick Ozark

"""
Inputs: (1) A normCounts file from htseq
        (2) A 13-column peak annotation file with the nearest TSS for each observed peak
"""

import sys, getopt, os
def main():
    input_anno = ""
    input_counts = ""
    try:
        opts,args = getopt.getopt(sys.argv[1:], "ha:c:",["help","annoFile=","countsFile="])
    except getopt.GetoptError:
        print("python3 chipPeakToGeneExpression.py -a|--annoFile <peak anno file> -c|--countsFile <normalized gene counts from htseq>")
        sys.exit()
    for opt,arg in opts:
        if opt in ("-h","--help"):
            print("python3 chipPeakToGeneExpression.py -a|--annoFile <peak anno file> -c|--countsFile <normalized gene counts from htseq>")
            sys.exit()
        elif opt in ("-a","--annoFile"):
            input_anno = arg
        elif opt in ("-c","--countsFile"):
            input_counts = arg
        else:
            assert False, "Unhandled option"

    # Dictionary to store gene ID's as keys and normalized counts as values
    d = {}

    cdt_output = str(input_anno + ".cdt")
    anno = open(input_anno,'r')
    counts = open(input_counts,'r')
    cdt = open(cdt_output,'w')

    # Stores value of number of samples in counts file
    length_counts = 0
    counts_line_count = 1

    # Stores the names of each sample from the counts file
    column_names = ""

    # Loop through counts file and store counts in dictionary with the key as gene ID and value as normalized counts
    for line in counts:
        if counts_line_count == 1:
            length_counts = len(line.split()) - 1 # Subtract 1 because the "Gene" field doesn't count as a sample
            counts_line_count += 1
            column_names = line.split()[1:]
            pass
        else:
            line = line.split()
            line[0] = line[0].strip('\"')
            d[line[0]] = line[1:]
    counts.close()

    # Add header to cdt file
    cdt.write("Gene_id" + "\t" + "Gene_symbol" + "\t")
    for field in column_names:
        field = field.strip('\"')
        cdt.write(str(field) + "\t")
    cdt.write("\n")

    # line_count makes sure the column names are not stored as data in the cdt output
    anno_line_count = 1
    for line in anno:
        if anno_line_count == 1:
            anno_line_count += 1
            pass
        else:
            _chr, start, end, name, score, distToNearestGene, nearestGene, nearestGeneName, nearestGeneBiotype, distToNearestTSS, nearestTSS, tssGeneName, tssGeneBiotype = line.split()

            # Remove quotes from input anno file

#            _chr = _chr.strip('\"')
#            name = name.strip('\"')
#            score = score.strip('\"')
            nearestGene = nearestGene.strip('\"')
            nearestGeneName = nearestGeneName.strip('\"')
#            nearestGeneBiotype = nearestGeneBiotype.strip('\"')
            nearestTSS = nearestTSS.strip('\"')
            tssGeneName = tssGeneName.strip('\"')
#            tssGeneBiotype = tssGeneBiotype.strip('\"')

            cdt.write(nearestTSS + "\t" + tssGeneName)

            if nearestTSS in d.keys():
                for field in d[nearestTSS]:
                    cdt.write("\t" + str(field))
                cdt.write("\n")
            else:
                for i in range(length_counts):
                    cdt.write("\t" + str(0))
                cdt.write("\n")

    anno.close()
    cdt.close()

if __name__ == "__main__":
    main()
