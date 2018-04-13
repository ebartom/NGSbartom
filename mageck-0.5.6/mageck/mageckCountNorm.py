""" MAGeCK count normalization related functions 
"""

from __future__ import print_function
import sys
import math
import logging
import string

def mageckcount_gettotalnormfactor(ctable):
  """
  Get the factor by total normalization
  """
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  # calculate the sum
  sumsample=[0]*n
  for (k,v) in ctable.items():
    sumsample=[sumsample[i]+v[i] for i in range(n)]
  # normalizing factor
  avgsample=sum(sumsample)/float(n)
  samplefactor=[avgsample/k for k in sumsample]
  return samplefactor

def mageckcount_getmediannormfactor(ctable):
  """
  Get the factor by median normalization
  """
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  meanval={k:math.exp( (sum( [ math.log(v2+1.0) for v2 in v])*1.0/n) ) for (k,v) in ctable.items() if sum(v)>0}  # geometric mean
  meanval={k:(lambda x: x if x>0 else 1)(v) for (k,v) in meanval.items()} # delete those with all 0 read counts
  #samplefactor=[0]*n
  medianfactor=[0.0]*n
  for ni in range(n):
    meanfactor=[ v[ni]/meanval[k] for (k,v) in ctable.items() if k in meanval]
    #print(str(sorted(meanfactor)))
    xfactor=sorted(meanfactor)[len(meanfactor)//2] # corrected 
    if xfactor>0.0:
      medianfactor[ni]=1.0/xfactor
      #logging.debug('xfactor:'+str(xfactor))
  return medianfactor



def normalizeCounts(ctable,sgdict=None,method='median',returnfactor=False,reversefactor=False,controlsgfile=None):
  """
  Central function for normalizing read counts
  Parameters:
  --------------
  ctable
    A dictionary of read counts: {sgrnaseq:[count0,count1,...]} 
  sgdict
    A dictionary of sgRNAs: {sgrnaseq: (sgrnaid, geneid)} (for count table input, sgrnaid=sgrnaseq)
  method
    Normalization methods: none,total,median,control
  returnfactor
    Whether to normalize read counts, or just return factors
  reversefactor
    Whether the factor should be reversed (1/factor)
  controlsgfile
    A file name containing control sgRNAs
 
  Return value: 
  --------------
  {sgRNA:[read counts]} if returnfactor == False, or [size_factor] if returnfactor == True
  
  By default, for higher read depths, the size factor is <1. If reversefactor is set, the factor is >1 (or 1/factor) 
  """
  # sums
  if len(ctable)==0:
    return ctable.copy()
  n=len(ctable[list(ctable.keys())[0]]) # samples
  m=len(ctable) # sgRNAs
  # calculate the total factor
  if method=='control':
    # use control sgRNAs 
    if controlsgfile == None:
      logging.error('Error: a list of control sgRNAs should be specified when using control normalization.')
      sys.exit(-1)
    controlsglist=[line.strip() for line in open(controlsgfile)]
    logging.info('Loaded '+str(len(controlsglist))+' control sgRNAs from '+controlsgfile)
    # get the tables for control sgRNAs
    ctable_nm={k:v for (k,v) in ctable.items() if k in controlsglist}
    if len(ctable_nm)==0 and sgdict !=None: # try to search in sgdict
      ctable_nm={k:v for (k,v) in ctable.items() if sgdict[k][0] in controlsglist}
    logging.info('control sgRNAs for normalization:'+str(len(ctable_nm)))
    if len(ctable_nm) ==0:
      logging.error('Error: cannot find records of control sgRNAs in the read count table.')
      sys.exit(-1) 
    method='median' # force to use median of controls
  else:
    ctable_nm=ctable
  samplefactor=mageckcount_gettotalnormfactor(ctable_nm)
  logging.debug('Initial (total) size factor: '+' '.join([str(x) for x in samplefactor]))
  if method=='median':
    # calculate the medianfactor
    medianfactor=mageckcount_getmediannormfactor(ctable_nm)
    usetotalnorm=False
    for ni in range(n):
      if medianfactor[ni]==0.0:
        logging.warning('Sample '+str(ni)+' has zero median count, so median normalization is not possible. Switch to total read count normalization.')
        usetotalnorm=True
    if usetotalnorm == False:
      samplefactor=medianfactor
      logging.debug('Median factor: '+' '.join([str(x) for x in samplefactor]))
    # end median normalization
  elif method=='none': 
    # no normalization
    samplefactor=[1.0]*n
  logging.info('Final size factor: '+' '.join([str(x) for x in samplefactor]))
  
  #
  if returnfactor:
    if reversefactor:
      return [1.0/x for x in samplefactor]
    else:
      return samplefactor
  # normalize the table
  ntable={ k: [ samplefactor[i]*v[i] for i in range(n)] for (k,v) in ctable.items()}
  return ntable
  


