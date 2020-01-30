#
#
from __future__ import print_function

import sys
import re
import os
import logging
from mageck.fileOps import *
from mageck.mageckCount import *


class VisualRValue:
  '''
  Class for R visualization
  '''
  outprefix='sample1'
  genesummaryfile=''
  cpindex=[]
  targetgene=[]
  cplabel=''

  # internal variable, for R file
  rtemplatestr='';  # template string
  rtemplate_gene_str='';  # template string
  outrfh=None;  # file handle for R file

  # for Rnw file
  rnwtemplatestr=''
  outrnwfh=None
  outrnwstring=''
  # for statistics of gene_summary_file
  comparisonlabel=[]; # label for comparison
  ngenes=[]; # number of genes
  selection=[]; # selections
  nfdr1=[]; # genes with FDR < 1, 5, 25%
  nfdr5=[]
  nfdr25=[]
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
    filename=os.path.join(os.path.dirname(__file__),'plot_template.RTemplate')
    if os.path.isfile(filename) and os.path.exists(filename):
      logging.info('Loading R template file: '+filename+'.')
    else:
      logging.error('Cannot find template file: '+filename)
      return -1
    filename_indgene=os.path.join(os.path.dirname(__file__),'plot_template_indvgene.RTemplate')
    if os.path.isfile(filename_indgene) and os.path.exists(filename_indgene):
      logging.info('Loading R template file: '+filename_indgene+'.')
    else:
      logging.error('Cannot find template file: '+filename_indgene)
      return -1
    # Rnw files
    filename_rnw=os.path.join(os.path.dirname(__file__),'plot_template.Rnw')
    if os.path.isfile(filename_rnw) and os.path.exists(filename_rnw):
      logging.info('Loading Rnw template file: '+filename_rnw+'.')
    else:
      logging.error('Cannot find template file: '+filename_rnw)
      return -1
    logging.debug('Setting up the visualization module...')
    #
    # loading
    with open(filename, "r") as rtfile:
      rtp=rtfile.read()
      outpdffile=self.outprefix+'.pdf'
      # rtp=re.sub('__OUTPUT_FILE__',outpdffile,rtp); # pdf file
      #
      outrfile=self.outprefix+'.R'

      # write to R file
      outrfh=open(outrfile,'w')
      self.outrfh=outrfh
      self.rtemplatestr=rtp

      # write pdf loading
      pdfloadstr="pdf(file='"+os.path.basename(outpdffile)+"',width=4.5,height=4.5);"
      # write file reading
      # rtp=re.sub('__GENE_SUMMARY_FILE__',self.genesummaryfile,rtp); # gene_summary
      tableadstr="gstable=read.table('"+os.path.basename(self.genesummaryfile)+"',header=T)"
      print(pdfloadstr,file=outrfh)
      print(tableadstr,file=outrfh)
    # load individual gene code
    with open(filename_indgene, "r") as rtfile:
      rtp=rtfile.read()
      self.rtemplate_gene_str=rtp
    # load Rnw file
    with open(filename_rnw,"r") as rtfile:
      rnw=rtfile.read()
      self.rnwtemplatestr=rnw
      outrfile=self.outprefix+'_summary.Rnw'
      self.outrnwstring=self.rnwtemplatestr
      outrfh=open(outrfile,'w')
      self.outrnwfh=outrfh

    return 0

  def closeRTemplate(self):
    '''
    Close the R file
    '''
    # write to R file
    print("dev.off()",file=self.outrfh)
    #
    rnwfile=self.outprefix+'_summary.Rnw'
    rfile=self.outprefix+'.R'
    summaryfile=self.outprefix+'_summary'
    latexfile=self.outprefix+'_summary.tex'
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile)
    # write code in R file to generate PDF files
    print("Sweave(\""+rnwfile_base+"\");\nlibrary(tools);\n",file=self.outrfh)
    print("texi2dvi(\""+os.path.basename(latexfile)+"\",pdf=TRUE);\n",file=self.outrfh)
    # write to Rnw file
    print(self.outrnwstring,file=self.outrnwfh)

    self.outrnwfh.close()
    self.outrfh.close()

  def WriteRTemplate(self):
    '''
    Given a VisualRValue object, write to an R file.
    The following variables need to be set up: self.cpindex, self.targetgene, self.cplabel
    '''

    # load file
    rtp=self.rtemplatestr
    # replace the variables


    indexchr=','.join([str(x) for x in self.cpindex])
    rtp=re.sub('__INDEX_COLUMN__',indexchr,rtp); # index
    logging.debug('Column index:'+indexchr)

    targetchr="c("+','.join(['"'+x+'"' for x in self.targetgene])+")"
    rtp=re.sub('__TARGET_GENE_LIST__',targetchr,rtp); # index

    rtp=re.sub('__SAMPLE_LABEL__',self.cplabel,rtp)

    # write to R file
    print(rtp,file=self.outrfh)
    # save to Rnw file
    rtprnw=""
    rtprnw+=r"\n\\newpage\\section{Comparison results of "+re.sub('_',' ',self.cplabel)+"}\n"
    rtprnw+=r"\n"+"The following figure shows the distribution of RRA score in the comparison "+re.sub('_',' ',self.cplabel)+", and the RRA scores of "+str(len(self.targetgene))+" genes.\n"
    rtprnw+="\n<<echo=FALSE>>=\n"
    tableadstr="gstable=read.table('"+os.path.basename(self.genesummaryfile)+"',header=T)"
    rtprnw+=tableadstr+"\n@"+"\n"
    rtprnw+=r"%\n\n\n"
    # rtprnw+=r"\\"+"begin{figure}[!h]\n"+r"\\"+"begin{center}\n"
    rtprnw+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>="+rtp+"@"+"\n"
    # rtprnw+="\\end{center}\n\\end{figure}\n"
    rtprnw+=r"%%\n"+"\\clearpage\n"
    rtprnw+="%__INDIVIDUAL_PAGE__\n\n"
    updatestr=re.sub('%__INDIVIDUAL_PAGE__',rtprnw,self.outrnwstring)
    self.outrnwstring=updatestr

  def loadTopK(self, filename, k=10):
    '''
    Load the top k gene names from the file
    '''
    n=-1
    self.targetgene=[]
    for line in open(filename):
      n+=1
      if n==0:
        continue
      if n<=k:
        field=line.strip().split()
        tgenename=field[0]
        self.targetgene+=[tgenename]
      else:
        break
    # write to file?
    logging.info('Loading top '+str(k) +' genes from '+filename+': '+','.join(self.targetgene))
    self.WriteRTemplate()
    return 0

  def loadTopKWithExp(self,filename,nttab,sgrna2genelist,collabels,k=10):
    '''
    Plot the individual sgRNA read counts of top k genes, and the position of these gene scores
    '''
    self.loadTopK(filename,k)
    self.loadGeneExp(self.targetgene,nttab,sgrna2genelist,collabels)

  def loadSelGeneWithExp(self,targetgene,nttab,sgrna2genelist,collabels,k=10):
    '''
    Plot the individual sgRNA read counts of top k genes, and the position of these gene scores
    '''
    self.targetgene=targetgene
    self.WriteRTemplate()
    self.loadGeneExp(self.targetgene,nttab,sgrna2genelist,collabels)

  def loadSelGene(self,targetgene):
    self.targetgene=targetgene
    self.WriteRTemplate()

  def loadGeneExp(self,genelist,nttab,sgrna2genelist,collabels):
    '''
    Load the sgRNA read counts of selected genes into file
    '''
    # insertion str
    insertstr=''
    # set up par
    #
    #parstr="<<echo=FALSE>>=\n"+" par(mfrow=c(2,2));\n" +"@"+"\n"
    #insertstr+=parstr
    nsubfigs=4
    npl=0
    # explanation
    insertstr+=r"\\newpage\n"+"The following figures show the distribution of sgRNA read counts (normalized) of selected genes in selected samples.\n"
    # plot individual genes
    for gene in genelist:
      sglist=[ k for (k,v) in sgrna2genelist.iteritems() if v==gene]
      ntgene={k:v for (k,v) in nttab.iteritems() if k in sglist}
      npl+=1
      # load to file
      valstring='list('
      vstrlist=[]
      for (k,v) in ntgene.iteritems():
        vstr='c('+','.join([str(vv) for vv in v])+')'
        vstrlist+=[vstr]
      valstring+=','.join(vstrlist)
      valstring+=')'
      rtp=self.rtemplate_gene_str
      rtp=re.sub('__TARGET_GENE__','"'+gene+'"',rtp)
      rtp=re.sub('__TARGET_MATRIX__',valstring,rtp)
      # labels
      clabstr='c('
      clabstr+=','.join(['"'+ x+'"' for x in collabels])
      clabstr+=')'
      rtp=re.sub('__COL_LABEL__',clabstr,rtp)
      # save to R file
      print(rtp,file=self.outrfh)
      # save to Rnw file
      rtprnw=''
      if npl %4 ==1:
        rtprnw=r"%\n\n\n"
        # rtprnw+=r"\\"+"begin{figure}[!h]\n"+r"\\"+"begin{center}\n"
        rtprnw+="<<fig=TRUE,echo=FALSE,width=4.5,height=4.5>>=\n"
        rtprnw+="par(mfrow=c(2,2));\n"
      rtprnw+=rtp
      if npl%4==0 or npl == len(genelist):
        rtprnw+="\npar(mfrow=c(1,1));\n"
        rtprnw+="@"+"\n"
        # rtprnw+="\\end{center}\n\\end{figure}\n%%\n"
      insertstr+=rtprnw
      # rtprnw+="%__INDIVIDUAL_PAGE__\n\n"
      # updatestr=re.sub('%__INDIVIDUAL_PAGE__',rtprnw,self.outrnwstring)
      # self.outrnwstring=updatestr
    # recover par
    ## parstr="<<echo=FALSE>>=\n"+" par(mfrow=c(1,1));\n" +"@"+"\n"
    ## insertstr+=parstr
    # write to Rnw file
    insertstr+="%__INDIVIDUAL_PAGE__\n\n"
    updatestr=re.sub('%__INDIVIDUAL_PAGE__',insertstr,self.outrnwstring)
    self.outrnwstring=updatestr

  def getGeneSummaryStat(self,args,isplot=True):
    '''
    Get the summary statistics of gene summary file
    '''
    n=0
    ncomparisons=0
    comparisonlabel=self.comparisonlabel
    for line in open(self.genesummaryfile):
      n+=1
      field=line.strip().split('\t')
      if n==1:
        if len(field) %12 !=2:
          logging.error('Not enough field in gene summary file: '+args.gene_summary)
          sys.exit(-1)
        ncomparisons=int( (len(field)-2)/12)
        # extract comparison labels
        for i in range(ncomparisons):
          neglabelindex=i*12+2
          negstr=re.sub('.lo.neg','',field[neglabelindex])
          comparisonlabel+=[negstr]
          comparisonlabel+=[negstr]
        # set up the variables
        self.nfdr1=[0]*2*ncomparisons
        self.nfdr5=[0]*2*ncomparisons
        self.nfdr25=[0]*2*ncomparisons
      else:
        for i in range(ncomparisons):
          nneg=i*12+2+2; npos=i*12+2+8
          try:
            if float(field[nneg])<0.01:
              self.nfdr1[2*i]+=1
            if float(field[npos])<0.01:
              self.nfdr1[2*i+1]+=1
            if float(field[nneg])<0.05:
              self.nfdr5[2*i]+=1
            if float(field[npos])<0.05:
              self.nfdr5[2*i+1]+=1
            if float(field[nneg])<0.25:
              self.nfdr25[2*i]+=1
            if float(field[npos])<0.25:
              self.nfdr25[2*i+1]+=1
          except ValueError:
            pass
      # end if
    # end for
    self.ngenes=[n-1]*2*ncomparisons
    self.selection=['negative','positive']*ncomparisons
    #
    if isplot==True:
      self.writeGeneSummaryStatToBuffer()

  def writeGeneSummaryStatToBuffer(self):
    '''
    Write statistics from gene summary file to buffer
    '''
    # insert string
    insertstr=''
    insertstr+='comparisons=c(' + ','.join(['"'+x+'"' for x in self.comparisonlabel ])  +');\n'
    insertstr+='ngenes=c('+ ','.join([str(x) for x in self.ngenes]) +');\n'
    insertstr+='direction=c('+','.join(['"'+x+'"' for x in self.selection])+');\n'
    insertstr+='fdr1=c('+','.join([str(x) for x in self.nfdr1])+');\n'
    insertstr+='fdr5=c('+','.join([str(x) for x in self.nfdr5])+');\n'
    insertstr+='fdr25=c('+','.join([str(x) for x in self.nfdr25])+');\n'
    #
    nwktowrite=re.sub('#__GENE_SUMMARY_STAT__',insertstr,self.outrnwstring)
    self.outrnwstring=nwktowrite

  def generatePDF(self,keeptmp=False):
    '''
    Call R and pdflatex
    '''
    rnwfile=self.outprefix+'_summary.Rnw'
    rfile=self.outprefix+'.R'
    summaryfile=self.outprefix+'_summary'
    (rnwfile_dir,rnwfile_base)=os.path.split(rnwfile)
    if rnwfile_dir=='':
      rnwfile_dir='./'
    systemcall('cd '+rnwfile_dir+'; '+'Rscript '+os.path.basename(rfile))
    #systemcall('cd '+rnwfile_dir+'; '+ 'R CMD Sweave '+rnwfile_base)
    #systemcall('export SWEAVE_STYLEPATH_DEFAULT="TRUE";'+ 'cd '+rnwfile_dir+'; '+'pdflatex '+os.path.basename(summaryfile))
    # cleaning the fraction pdf
    if keeptmp==False:
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'-*.pdf')
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'.aux')
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'.tex')
      systemcall('cd '+rnwfile_dir+'; '+'rm -rf '+os.path.basename(summaryfile)+'.toc')





def plot_main(args):
  '''
  Main entry for plotting
  '''
  # loading count tables
  mapres=getcounttablefromfile(args.count_table)
  cttab=mapres[0]
  sgrna2genelist=mapres[1]
  samplelabelindex=mapres[2]

  # parse labels
  (treatgroup,treatgrouplabellist)=parse_sampleids(args.samples,samplelabelindex)
  # parse selected genes
  if args.genes==None:
    selgene=[]
  else:
    selgene=args.genes.split(',')

  # initialize R visualization init
  vrv=VisualRValue()
  vrv.setPrefix(args.output_prefix)
  vrv.genesummaryfile=args.gene_summary
  vrv.startRTemplate()

  # generate summary file; must be done before plotting any individual genes
  vrv.getGeneSummaryStat(args)

  # check the maximum column in gene summary
  n=0
  ncomparisons=0
  comparisonlabel=[]
  for line in open(vrv.genesummaryfile):
    n+=1
    if n==1:
      field=line.strip().split('\t')
      if len(field) %12 !=2:
        logging.error('Not enough field in gene summary file: '+args.gene_summary)
        sys.exit(-1)
      ncomparisons=int( (len(field)-2)/12)
      # extract comparison labels
      for i in range(ncomparisons):
        neglabelindex=i*12+2
        negstr=re.sub('.lo.neg','',field[neglabelindex])
        comparisonlabel+=[negstr]
    else:
      break

  # read the sgRNA-gene table for rank association
  # normalization
  cttab_sel={k:([v[i] for i in treatgroup]) for (k,v) in cttab.iteritems()}; # controlgroup do not overlap with treatgroup
  if hasattr(args,'norm_method'):
    nttab=normalizeCounts(cttab_sel,method=args.norm_method,controlsgfile=args.control_sgrna)
  else:
    nttab=normalizeCounts(cttab_sel)

  if len(selgene)>0:
    vrv.loadGeneExp(selgene,nttab,sgrna2genelist,treatgrouplabellist)
  # testing the comparisons
  for nc in range(ncomparisons):
    # visualization: load top k genes
    # print(str(samplelabelindex))
    vrv.cplabel=comparisonlabel[nc]+' neg.'
    vrv.cpindex=[2+12*nc+1]
    vrv.loadSelGene(selgene)
    vrv.cplabel=comparisonlabel[nc]+' pos.'
    vrv.cpindex=[2+12*nc+6+1]
    vrv.loadSelGene(selgene)




  # generate pdf file
  vrv.closeRTemplate()
  vrv.generatePDF(args.keep_tmp)
  #systemcall('Rscript '+vrv.outprefix+'.R')
  #systemcall('R CMD Sweave '+vrv.outprefix+'_summary.Rnw')
  #systemcall('pdflatex '+vrv.outprefix+'_summary')
