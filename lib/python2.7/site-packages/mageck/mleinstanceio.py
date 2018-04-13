'''
Reading and writing instances
'''

from __future__ import print_function
import sys
import numpy as np
from mageck.mleclassdef import *

import logging

def read_gene_from_file(filename,includesamples=None):
  '''
  Reading gene models 
  Parameters
  ----------
  filename
    file name of the read count table 
  includesamples
    If not None, only samples in the includesampels are included
  '''
  # first, read count table
  allgenedict={}
  nline=0
  nsamples=0
  ngene=0
  
  sampleids=[]
  sampleindex=[]
  sampleids_toindex={}
  hascsv=False
  if filename.upper().endswith('.CSV'):
    hascsv=True
    logging.info('Treating '+filename+' as csv file ...')
  for line in open(filename):
    nline+=1
    if hascsv==False:
      field=line.strip().split()
    else:
      field=line.strip().split(',')
    if nline==1:
      # The first line: check sample columns
      nsamples=len(field)-2
      sampleids=field[2:]
      for i in range(nsamples):
        sampleids_toindex[sampleids[i]]=i
      if includesamples != None:
        logging.info('Loaded samples:'+';'.join(includesamples))
        for si in includesamples:
          if si not in sampleids_toindex:
            logging.error('Sample '+si+' cannot be found on the original read count table '+filename)
            sys.exit(-1)
        sampleindex=[sampleids_toindex[si] for si in includesamples]
        logging.info('Sample index: '+';'.join([str(x) for x in sampleindex]))
      else:
        sampleindex=[i for i in range(nsamples)]
      continue
    sgid=field[0]
    gid=field[1]
    if gid not in allgenedict:
      sks=SimCaseSimple()
      sks.prefix=gid
      sks.nb_count=[]
      sks.sgrnaid=[]
      ngene+=1
      for i in sampleindex:
        sks.nb_count+=[[]]
      allgenedict[gid]=sks
    else:
      sks=allgenedict[gid]
    sks.sgrnaid+=[sgid]
    for i in range(len(sampleindex)): 
      ni=sampleindex[i]
      try:
        nrt=float(field[ni+2])+1 # add 1 pseudocount
        sks.nb_count[i]+=[nrt]
      except ValueError:
        print('Error loading line '+str(nline))
  # end for loop
  logging.info('Loaded '+str(ngene)+' genes.')
  #
  # convert nb_count to matrix
  for (gene,ginst) in allgenedict.iteritems():
    ginst.nb_count=np.matrix(ginst.nb_count)
  return allgenedict

def write_gene_to_file(allgenedict,outfile,args,betalabels=None):
  '''
  Write gene to file
  '''
  ofid=open(outfile,'w')
  tmpinst=allgenedict[allgenedict.keys()[0]]
  nbeta=len(tmpinst.beta_estimate)-(tmpinst.nb_count.shape[1])
  # print header
  # headerterms=['|beta','|z','|neg|p-value','|neg|fdr','|pos|p-value','|pos|fdr', '|neg|permutation', '|neg|permutation-fdr', '|pos|permutation','|pos|permutation-fdr' ] # one-sided
  # headerterms=['|beta','|z','|p-value','|fdr', '|permutation', '|permutation-fdr' ] # two-sided,using Wald test for p value
  headerterms=['|beta','|z','|p-value','|fdr','|wald-p-value','|wald-fdr' ] # two-sided,using permutation test for p value
  if betalabels == None: 
    # reportlabels='\t'.join(['\t'.join(['beta_'+str(i+1)+'|beta','beta_'+str(i+1)+'|z','beta_'+str(i+1)+'|p-value','beta_'+str(i+1)+'|permutation']) for i in range(nbeta)])  # two-sided
    reportlabels='\t'.join(['\t'.join(['beta_'+str(i+1)+ hst for hst in headerterms ]) for i in range(nbeta)]) # one-sided
  else:
    if len(betalabels)-1!=nbeta:
      raise ValueError('Beta labels do not match to columns of the design matrix.')
    # reportlabels='\t'.join(['\t'.join([sstt+'|beta',sstt+'|z',sstt+'|p-value',sstt+'|permutation']) for sstt in betalabels[1:]]) # two-sided
    reportlabels='\t'.join(['\t'.join([sstt+hst for hst in headerterms]) for sstt in betalabels[1:]])
  print('\t'.join(['Gene','sgRNA',reportlabels]),file=ofid)
  # print for each gene 
  if args.control_sgrna is not None:
    controlsglist=[line.strip() for line in open(args.control_sgrna)]
  else:
    controlsglist=[]
  for (tgid,tginst) in allgenedict.iteritems():
    # skip genes if it consists of all control sgRNAs
    notskip=False
    for sgid in tginst.sgrnaid:
      if sgid not in controlsglist:
        notskip=True
    if notskip==False:
      continue
    wfield=tginst.gene_to_printfield()
    print('\t'.join(wfield),file=ofid)
  # end
  ofid.close()

def write_sgrna_to_file(allgenedict,outfile):
  '''
  Write gene to file
  '''
  ofid=open(outfile,'w')
  tmpinst=allgenedict[allgenedict.keys()[0]]
  nbeta=len(tmpinst.beta_estimate)-len(tmpinst.nb_count[0])
  print('\t'.join(['Gene','sgRNA','eff']),file=ofid)
  for (tgid,tginst) in allgenedict.iteritems():
    for i in range(len(tginst.w_estimate)):
      wfield=[tginst.prefix,tginst.sgrnaid[i],decformat(tginst.w_estimate[i])]
      print('\t'.join(wfield),file=ofid)
  ofid.close()





