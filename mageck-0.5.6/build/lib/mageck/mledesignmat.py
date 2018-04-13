'''
Defining the design matrix
'''
from __future__ import print_function

import numpy as np
import logging
import sys

def getsimpledesignmat(n):
  '''
  For a two smaple case, get a design matrix
  '''
  x1=np.identity(n)
  x2=np.zeros(n)
  x3=np.ones(n)
  c1=np.column_stack((x1,x2))
  c2=np.column_stack((x1,x3))
  c3=np.vstack((c1,c2))
  return np.matrix(c3)

def getextenddesignmat(nsg,nsample,designmat,includebase=False):
  '''
  Get an extended design matrix from an arbitrily design matrix
  Parameters:
    nsg
        number of sgRNAs
    nsample
        number of samples (excluding baseline samples)
    includebase
        whether to include further nsgrna*nsample rows at the end
  Return value:
    extended design matrix

  the matrix size:
    (nsg* (nsample+1))*(nbeta) if includebase = False
    (nsg* (2*nsample+1))*(nbeta) if includebase = True
  '''
  if designmat.shape[0] !=nsample:
    raise ValueError('The row of design matrix does match the number of samples.')

  sgiden=np.identity(nsg)
  matright=np.matrix([[0]*designmat.shape[1]]*nsg)
  c3_top=np.column_stack((sgiden,matright))
  c3=c3_top.copy()
  for ns in range(nsample):
    # expand the ns-th row of design mat
    matright=np.matrix(designmat[ns].tolist()*nsg)
    px=np.column_stack((sgiden,matright))
    c3=np.vstack((c3,px))
  if includebase:
    for ns in range(nsample):
      c3=np.vstack((c3,c3_top))
  return c3


def analyze_designmat(design_mat):
  '''
  analyze the design matrix
  Return value:
    the index of baseline samples
    the new design matrix
  '''
  init_design_mat=design_mat
  init_design_mat_rowsum=[x[0] for x in np.sum(init_design_mat,axis=1).tolist()]
  for i in range(len(init_design_mat_rowsum)):
    if init_design_mat_rowsum[i]==0:
      print('Error parsing the design matrix: 0 row sums for some samples.')
  init_design_basesample=[t for t in range(init_design_mat.shape[0]) if init_design_mat_rowsum[t]==1]
  init_design_othersample=[t for t in range(init_design_mat.shape[0]) if init_design_mat_rowsum[t]>1]
  new_design_matrix=design_mat[init_design_othersample,1:]
  return (init_design_basesample,new_design_matrix)

def parse_designmat_from_day0(args):
  '''
  Parse the design matrix from a string.
  Parameters
  ----------
  args 
    An argparse object where day0_label and count_table will be used

  Return values
  ----------
  ndesm
    The design matrix object
  samples
    Sample id (rows of design matrix)
  beta_label
    Beta label (columns of design matrix)
  '''
  # first, get the labels in the count table
  if args.count_table.upper().endswith('.CSV'):
    splitter=','
  else:
    splitter=None
  with open(args.count_table) as f:
    firstline=f.readline()
  field=firstline.strip().split(splitter)
  samples=field[2:]

  # check day0 label
  args.day0_label=args.day0_label.split(',')
  for dl in args.day0_label:
    if dl not in samples:
      logging.error('Label '+dl+' specified in the --day0-label option does not match count table. Please double check.')
      sys.exit(-1)
  # check day0s
  nsample=len(samples)
  nonday0sample=[x for x in samples if x not in args.day0_label]
  n_nond0=len(nonday0sample)
  if n_nond0 == nsample:
    logging.error('Cannot find day0 sample.')
    sys.exit(-1)
  sample_order=args.day0_label+nonday0sample
  # dessplit=[ [1 ] + [(lambda t:1 if t in nonday0sample else 0)(tlb) for tlb in sample_order]]
  dessplit=[ [1]+ [  (nd==tlb)*1 for nd in nonday0sample ]   for tlb in sample_order ]
  desmat=np.matrix(dessplit)
  return (desmat,sample_order,['baseline']+nonday0sample)
    


def parse_designmat(designmatstr):
  '''
  Parse the design matrix from a string.
  Parameters
  ----------
  designmatstr
    A string of design matrix, either specifying the matrix elements directly, or a file containing the design matrix. See explanations below.

  Return values
  ----------
  ndesm
    The design matrix object
  samples
    Sample id (rows of design matrix)
  beta_label
    Beta label (columns of design matrix)
  ----------
  Two different values are specified: a file name containing the design matrix, or a string of design matrix.

  If designmatstr is a string, use ";" to separate rows and "," to separate columns.
  For example, if a design matrix is as follows:

  1 0 0
  1 0 0
  1 1 0
  1 0 1

  Then the design matrix string should be:

  1,0,0;1,0,0;1,1,0;1,0,1

  Otherwise if a file is specified, the file should be a tab-delimited file containing the design matrix. Row and column names are required:
  sample b0 b1 b2
  s1     1  0  0
  s2     1  0  0
  s3     1  1  0
  s4     1  0  1
  '''
  #
  (nret,ndesm)=parse_designmat_str(designmatstr)
  if nret==0:
    return (ndesm,None,None)
  else:
    logging.info('Cannot parse design matrix as a string; try to parse it as a file name ...')
    (nret2,ndesm2,rows,columns)=parse_designmat_file(designmatstr)
    if nret2==0:
      return (ndesm2,rows,columns)
    else:
      raise ValueError('parsing the design matrix from file '+designmatstr)
  return (None,None,None)


def parse_designmat_str(designmatstr):
  '''
  Parse the design matrix from a string
  For explanations, see parse_designmat()
  Return value:
    (x, design_mat)
    x:
      0 if success, -1 if not
    designmat:
      design matrix if success, None if not
  '''
  desfield=designmatstr.split(';')
  try:
    dessplit=[ [float(x) for x in fx.split(',')] for fx in desfield]
    dessplitmat=np.matrix(dessplit)
    return (0, dessplitmat)
  except ValueError:
    return (-1, None)

def parse_designmat_file(designmatstr):
  '''
  Parse the design matrix from a file
  For explanations, see parse_designmat()
  Return value:
  ------------
    (x, design_mat,rows,columns)
    x:
      0 if success, -1 if not
    designmat:
      design matrix if success, None if not
    rows:
      labels for rows (i.e., sample labels)
    columns:
      labels for columns (i.e., sample columns)
  '''
  filefield=[x.strip().split() for x in open(designmatstr,'r')]
  try:
    dessplit=[ [float(x) for x in fx[1:]] for fx in filefield[1:]]
    dessplitmat=np.matrix(dessplit)
    collabels=filefield[0][1:]
    rowlabels=[ fx[0] for fx in filefield[1:] if len(fx)>0]
    return (0, dessplitmat,rowlabels,collabels)
  except ValueError:
    return (-1, None,None,None)


class DesignMatCache:
  '''
  Store the design matrix and relevant information of different # sgRNAs, so we don't have to construct it every time
  '''
  cache={}
  @staticmethod
  def save_record(orig_design_mat,n):
    '''
    Save the design mat of sgRNA n to cache
    '''
    # calculate the number of base line samples
    #
    nsample=orig_design_mat.shape[0]-1 # the number of samples excluding the 1st base sample
    (basesampleid,_new_designmat)=analyze_designmat(orig_design_mat)
    design_mat=orig_design_mat[1:,1:] # assuming the 1st column must be the baseline condition, and the 1st row is the base line condition; they are removed.
    extdesign_mat=getextenddesignmat(n,nsample,design_mat,includebase=True)
    extdesignmat_residule=extdesign_mat[0:(nsample+1)*n,n:] # this is the design matrix containing only non-baseline betas and first (nsample+1) samples
    # save
    save_tuple=(basesampleid,design_mat,extdesign_mat,extdesignmat_residule)
    DesignMatCache.cache[n]=save_tuple

  @staticmethod
  def has_record(n):
    '''
    If the cache has the record of design matrix of sgRNA n
    '''
    if n in DesignMatCache.cache:
      return True
    else:
      return False
  @staticmethod
  def get_record(n):
    '''
    Get the record in design matrix
    '''
    return DesignMatCache.cache[n]
