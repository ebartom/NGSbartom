# Author: Patrick Ozark
# Creates a script to plot fingerPrint of ChIPseq samples

import sys, getopt, os
def main():
    bam_path = ""
    out_path = ""
    script = "plot_fingerprint.sh"
    out = "fingerprint.pdf"
    try:
        opts,args = getopt.getopt(sys.argv[1:], "hi:o:",["help","input=","output="])
    except getopt.GetoptError:
        print("python3 getPlotHeatmap.py -i|--input <input BAM directory> -o|--output <output plot directory>")
        sys.exit()
    for opt,arg in opts:
        if opt in ("-h","--help"):
            print("python3 getPlotHeatmap.py -i|--input <input BAM directory> -o|--output <output plot directory>")
            sys.exit()
        elif opt in ("-i","--input"):
            bam_path = arg
        elif opt in ("-o","--output"):
            out_path = arg
        else:
            assert False, "Unhandled option"

    script_path = str(out_path) + "/" + str(script)
    output = open(script_path,"w")

    output.write("#!/bin/bash\n")
    output.write("#MSUB -A b1042\n")
    output.write("#MSUB -q genomics\n")
    output.write("#MSUB -l walltime=24:00:00\n")
    output.write("#MSUB -m a\n")
    output.write("#MSUB -j oe\n")
    output.write("#MOAB -W umask=0113\n")
    output.write("#MSUB -N ChIPfingerprint\n")
    output.write("#MSUB -l nodes=1:ppn=4\n")
    output.write("module load python/anaconda3.6\n")
    output.write("\n")
    output.write("plotFingerprint -b ")
    for file in os.listdir(bam_path):
        if file.endswith(".bam"):
            output.write(str(bam_path) + "/" + str(file) + " ")
    output.write("--smartLabels ")
    output.write("-p 4 ")
    output.write("--skipZeros ")
    output.write("--minMappingQuality 20 ")
    output.write("--plotTitle \"Fingerprint\" ")
    output.write("-o " + str(out_path) + "/" + str(out) + "\n")
    output.close()

if __name__ == "__main__":
    main()
