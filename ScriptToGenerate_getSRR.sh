#!/bin/sh
# Author: Pankhuri Wanjari
# Description: Create shell scripts to download SRR files using fastq-dump utility from SRA toolkit

echo "#!/bin/sh" > "getSRR.sh"
echo "#!/bin/sh \n\n# Export sra tool kit path" > "fastqdump.sh"
echo "export PATH=\$PATH:/projects/b1025/tools/sratoolkit.2.4.4-centos_linux64/bin" >> "fastqdump.sh"
echo >> "fastqdump.sh"
echo "#STEP1: Download SRR files using wget. The use fastq-dump from SRA toolkit to convert the SRA files to fastq format" >> "getSRR.sh"
echo "#STEP2: Concatenate fastq files from same samples and then rename it to an appropriate filename" >> "getSRR.sh"
echo "#STEP3: Gunzip the file" >> "getSRR.sh"
echo >> "getSRR.sh"
echo "#NOTE:\n#1. & is used to run the resp commands/processes in background\n#2. Wait is used to pause until execution of the background process has ended" >> "getSRR.sh"
echo >> "getSRR.sh"
#use `wc -l` command to get the number of lines/enteries in your input file
a=1

echo "# Script to convert SRR files to Fastq using fastq dump" >> "fastqdump.sh"

# get length of file
filelen=`cat SRRinfo.txt|wc -l`

while [ $a -le $filelen ]
do
   s=`expr $a`p
   cmd=`cat SRRinfo.txt | awk '{print $2}'| sed -n $s`
   cmd2=`cat SRRinfo.txt | awk '{print $2 ".fastq" " "  $3 ".fastq"}'| sed -n $s`
   cmd3=`cat SRRinfo.txt | awk '{print $3 ".fastq"}'| sed -n $s`
   echo "wget -P rawfiles ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP060/SRP060227/$cmd/$cmd.sra &" >> "getSRR.sh"
   echo "fastq-dump -O /projects/b1025/pwj/GSE70408/SRRFiles/rawfiles/fastqfiles/ rawfiles/$cmd.sra &" >> "fastqdump.sh"
   #echo "mv $cmd2" >> "getSRR.sh"
   #echo "gzip $cmd3 &" >> "getSRR.sh"
   r=`expr $a % 10`
   if [ $r == 0 ]
   then
      echo "wait" >> "getSRR.sh"
      echo >> "getSRR.sh"
      echo "wait" >> "fastqdump.sh"
      echo >> "fastqdump.sh"
   fi
   a=`expr $a + 1`
done

# *** ---- *** --- *** --- *** --- *** --- *** ---- *** --- *** --- *** --- *** --- *** ---- *** --- *** --- *** --- *** ---

# get unique sample names and store it in a text file
cat SRRinfo.txt | awk '{print $1}' | uniq > "temp.txt"
samplenum=`cat temp.txt|wc -l`
echo "$samplenum"

i=1
echo "#!/bin/sh" > "Rename.sh"
echo "#!/bin/sh" > "Gunzip.sh"


while [ $i -le $samplenum ]
do 
    s=`expr $i`p
    samplename=`cat temp.txt | awk '{print $1}'| sed -n $s`
    SRRpersample=`grep $samplename SRRinfo.txt | awk '{print $2}' | wc -l` #gives number of SRR files per sample
    j=1
    temp=()
    while [ $j -le $SRRpersample ] 
    do
        t=`expr $j`p
        cmdSRR=`grep $samplename SRRinfo.txt | awk '{print $2}' | sed -n $t`
        filename=`grep $samplename SRRinfo.txt | awk '{print $3}' | sed -n $t`
        temp[j]=$cmdSRR.fastq
        j=`expr $j + 1`
    done
    echo "cat ${temp[@]} > $filename.fastq" >> "Rename.sh"
    echo >> "Rename.sh"
    echo "gzip $filename.fastq &" >> "Gunzip.sh"
    r=`expr $i % 10`
    if [ $r == 0 ]
    then
       echo "wait" >> "Gunzip.sh"
       echo >> "Gunzip.sh"
    fi
    i=`expr $i + 1`
done


# Get number of SRR files for each sample and store in a variable
#set b=`grep "GSM1727070" SRRinfo.txt|awk '{print $2}'|wc -l`
#$b

#cat SRRinfo.txt|awk '{print $1}'|sed -n 1p
#grep "GSM1727070" SRRinfo.txt|awk '{print $2}'


# Two ways to download data

#wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP060/SRP060227/SRR2084433/SRR2084433.sra
#fastq-dump SRR2084433.sra

#srapath SRR084433
#you will get an output like this
#http://sra-download.ncbi.nlm.nih.gov/srapub/SRR2084431
#then use wget 
#wget -P Data http://sra-download.ncbi.nlm.nih.gov/srapub/SRR2084431