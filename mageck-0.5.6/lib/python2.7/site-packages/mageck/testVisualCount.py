#
#
from __future__ import print_function

import sys
import re
import os
import logging
from mageck.fileOps import *
from mageck.mageckCount import *


class VisualRCount:
  '''
  Class for generating reports of count command
  '''
  outprefix='sample1'

  # internal variable, for R file
  outrfh=None;  # file handle for R file

  # for Rnw file
  rnwtemplatestr=''
  outrnwfh=None
  outrnwstring=''
  # for statistics of coutns
  fastqfile=[] # fastq files
  fastqlabels=[] # fastq labels
  reads=[] # read counts
  mappedreads=[] # # mapped reads
  totalsgrnas=[] # # sgRNAs (in the library)
  zerocounts=[] # # 0-count sgRNAs
  gini=[]  # gini index
  # for QC
  negselqc=[] # measure the negative selection QC
  negselqcpval=[] # the p value
  negselqcpvalpermutation=[] # the p value by permutation
  negselqcpvalpermutationfdr=[] # the fdr value by permutation
  negselqcgene=[] # the number of genes for negative selection
  '''
  Member functions
  '''
  def setPrefix(self,prefix):
    '''
    Set up proper prefix
    '''
    (file_dir,file_base)=os.path.split(prefix)
    file_base_dot=re.sub('\.','_',file_base)
    self.outprefix=os.path.join(file_dir,file_base_dot)
  def startRTemplate(self):
    '''
    Open a template, create an R file

    '''
    # R files
    # Rnw files
    filename_rnw=os.path.join(os.path.dirname(__file__),'fastq_template.Rnw')
    if os.path.isfile(filename_rnw) and os.path.exists(filename_rnw):
      logging.info('Loading Rnw template file: '+filename_rnw+'.')
    else:
      logging.error('Cannot find template file: '+filename_rnw)
      return -1
    logging.debug('Setting up the visualization module...')
    #

    # R file
    outrfile=self.outprefix+'_countsummary.R'
    outrfh=open(outrfile,'w')
    self.outrfh=outrfh
    # load Rnw file
    with open(filename_rnw,"r") as rtfile:
      rnw=rtfile.read()
      self.rnwtemplatestr=rnw
      outrfile=self.outprefix+'_countsummary.Rnw'
      self.outrnwstring=self.rnwtemplatestr
      outrfh=open(outrfile,'w')
      self.outrnwfh=outrfh

    return 0

  def closeRTemplate(self):
    '''
    Close the R file
    '''
    # write to R file
    #
    rnwfile=self.outprefix+'_countsummary.Rnw'
    rfile=self.outprefix+'_countsummary.R'
    summaryfile=self.outprefix+'_countsummary'
    latexfile=self.outprefix+'_countsummary.tex'
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile)
    # write code in R file to generate PDF files
    print("Sweave(\""+rnwfile_base+"\");\nlibrary(tools);\n",file=self.outrfh)
    print("texi2dvi(\""+os.path.basename(latexfile)+"\",pdf=TRUE);\n",file=self.outrfh)
    # write to Rnw file
    print(self.outrnwstring,file=self.outrnwfh)

    self.outrnwfh.close()
    self.outrfh.close()




  def writeCountSummary(self):
    '''
    Write statistics from gene summary file to buffer
    '''
    # insert string
    insertstr=''
    import textwrap
    fastqwrap=[' '.join(textwrap.wrap(x,50)) for x in self.fastqfile]
    insertstr+='filelist=c(' + ','.join(['"'+x+'"' for x in fastqwrap])  +');\n'
    insertstr+='labellist=c('+ ','.join(['"'+x+'"' for x in self.fastqlabels]) +');\n'
    insertstr+='reads=c('+','.join([str(x) for x in self.reads])+');\n'
    insertstr+='mappedreads=c('+','.join([str(x) for x in self.mappedreads])+');\n'
    insertstr+='totalsgrnas=c('+','.join([str(x) for x in self.totalsgrnas])+');\n'
    insertstr+='zerocounts=c('+','.join([str(x) for x in self.zerocounts])+');\n'
    insertstr+='giniindex=c('+','.join([str(x) for x in self.gini])+');\n'
    #
    nwktowrite=re.sub('#__COUNT_SUMMARY_STAT__',insertstr,self.outrnwstring)
    # file names as list, instead of tables; disabled currently
    insertstr=''
    insertstr+='The fastq files processed are listed as follows.\n\n'
    insertstr+=r"\\"+"begin{enumerate}\n"
    for fq in self.fastqfile:
      insertstr+=r"\\"+"item "+fq+"\n"
    insertstr+=r"\\"+"end{enumerate}\n"
    # nwktowrite=re.sub('%__FILE_SUMMARY__',insertstr,nwktowrite)
    self.outrnwstring=nwktowrite

  def writeCountSummaryToTxt(self,txtfile):
    '''
    A stand-alone function to write the count summary to txt file
    fastqfile, fastqlabels, reads, mappedreads, and zerocounts must be set up
    '''
    ofstr=open(txtfile,'w')
    nsp=len(self.fastqfile)
    print('\t'.join(['File','Label','Reads','Mapped', 'Percentage','TotalsgRNAs',  'Zerocounts','GiniIndex','NegSelQC','NegSelQCPval', 'NegSelQCPvalPermutation','NegSelQCPvalPermutationFDR', 'NegSelQCGene']),file=ofstr)
    for i in range(nsp):
      print('\t'.join([self.fastqfile[i],
          self.fastqlabels[i],
          str(self.reads[i]),
          str(self.mappedreads[i]), 
          "{:.4g}".format(self.mappedreads[i]*1.0/self.reads[i]),
          str(self.totalsgrnas[i]),
          str(self.zerocounts[i]),
          "{:.4g}".format(self.gini[i]), # gini
          "{:.5g}".format(self.negselqc[i]), # QC terms
          "{:.5g}".format(self.negselqcpval[i]), # pval
          "{:.5g}".format(self.negselqcpvalpermutation[i]), # pval by permutation
          "{:.5g}".format(self.negselqcpvalpermutationfdr[i]), # pval by permutation
          str(self.negselqcgene[i])]),
          file=ofstr)
    ofstr.close()

  def insertReadCountBoxPlot(self,filename):
    '''
    Insert box plot of normalized read counts
    '''
    # inserted R code
    rtp="\ngenboxplot(\""+filename+"\");\n"
    #
    insertstr=''
    insertstr+=r"\n\\newpage\\section{Normalized read count distribution of all samples}\n"
    insertstr+="The following figure shows the distribution of median-normalized read counts in all samples.\n\n\n"
    insertstr+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n"
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    #
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring)
    #
    rtp="\ngenhistplot(\""+filename+"\");\n"
    insertstr='' 
    insertstr+="The following figure shows the histogram of median-normalized read counts in all samples.\n\n\n"
    insertstr+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n"
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,nwktowrite)
    #
    self.outrnwstring=nwktowrite
  #
  def insertClusteringPlot(self,filename):
    '''
    Insert clustering plot of normalized read counts
    '''
    # inserted R code
    rtp="\ngenclustering(\""+filename+"\");\n"
    #
    insertstr=''

    insertstr+=r"\n\\newpage\\section{Sample clustering}\n"
    insertstr+="The following figure shows the sample clustering result.\n\n\n"

    insertstr+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n"
    #
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring)
    self.outrnwstring=nwktowrite

  def insertPCAPlot(self,filename):
    '''
    Insert box plot of PCA analysis
    '''
    # inserted R code
    rtp="\ngenpcaplot(\""+filename+"\");\n"
    rtp2="\ngenpcavar("+");\n"
    #
    insertstr=''

    insertstr+=r"\n\\newpage\\section{Principle Component Analysis}\n"
    insertstr+="The following figure shows the first 2 principle components (PCs) from the Principle Component Analysis (PCA), and the percentage of variances explained by the top PCs.\n\n\n"

    insertstr+="\n<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n"
    insertstr+="\n<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp2+"@"+"\n"
    #
    insertstr+='\n%__INDIVIDUAL_PAGE__\n'
    nwktowrite=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring)
    self.outrnwstring=nwktowrite



  def generatePDF(self,keeptmp=False):
    '''
    Call R and pdflatex
    '''
    rnwfile=self.outprefix+'_countsummary.Rnw'
    rfile=self.outprefix+'_countsummary.R'
    summaryfile=self.outprefix+'_countsummary'
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile)
    if rnwfile_dir=='':
      rnwfile_dir='./'
    systemcall('cd '+rnwfile_dir+'; '+'Rscript '+os.path.basename(rfile))
    #systemcall('cd '+rnwfile_dir+'; '+ 'R CMD Sweave '+rnwfile_base)
    #systemcall('export SWEAVE_STYLEPATH_DEFAULT="TRUE";'+ 'cd '+rnwfile_dir+'; '+'pdflatex '+os.path.basename(summaryfile))
    # cleaning the fraction pdf
    if keeptmp==False:
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'-*.pdf')
