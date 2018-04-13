#!/usr/bin/env python
'''
MAGeCK MLE main entry

'''

from __future__ import print_function

import re
import sys
import random
import math
import logging

# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

# importing predefined functions
from mageck.mleargparse import mageckmle_postargs


def mageckmle_main(pvargs=None,parsedargs=None,returndict=False):
  '''
  Main entry for MAGeCK MLE
  ----
  Parameters:

  pvargs
    Arguments for parsing
  returndict
    If set true, will not try to run the whole prediction process, but will return after mean variance modeling
  '''
  # parsing arguments
  if parsedargs is not None:
    args=parsedargs
  else:
    args=mageckmle_parseargs(pvargs)
  args=mageckmle_postargs(args)
  # Bayes module
  if hasattr(args,'bayes') and args.bayes:
    from mlemageck_bayes import mageck_bayes_main
    sys.exit(0) # comment this when you think bayes module is completed
    mageck_bayes_main(parsedargs=args)
    sys.exit(0) # 
  # from mleclassdef import *
  # from mledesignmat import *
  # from mleem import *
  # from mleinstanceio import *
  # from mlemeanvar import *
  import scipy
  from scipy.stats import nbinom
  import numpy as np
  import numpy.linalg as linalg
  from mageck.mleinstanceio import read_gene_from_file,write_gene_to_file,write_sgrna_to_file
  from mageck.mleem import iteratenbem
  from mageck.mlemeanvar import MeanVarModel
  from mageck.mageckCount import normalizeCounts
  from mageck.mlesgeff import read_sgrna_eff,sgrna_eff_initial_guess
  from dispersion_characterization import sgrna_wide_dispersion_estimation_MAP_v2
  from mageck.mlemultiprocessing import runem_multiproc,iteratenbem_permutation
  from mageck.cnv_normalization import read_CNVdata,match_sgrnaCN,betascore_piecewisenorm,betascore_piecewisenorm

  # main process
  maxfittinggene=args.genes_varmodeling
  maxgene=np.inf
  # reading sgRNA efficiency
  read_sgrna_eff(args)
  # reading read count table
  allgenedict=read_gene_from_file(args.count_table,includesamples=args.include_samples)
  # 
  sgrna_eff_initial_guess(args,allgenedict)
  # calculate the size factor
  cttab_sel={}
  for (geneid,gk) in allgenedict.iteritems():
    sgid=gk.sgrnaid
    sgreadmat=gk.nb_count.getT().tolist()
    for i in range(len(sgid)):
      cttab_sel[sgid[i]]=sgreadmat[i]
  if hasattr(args,'norm_method'):
    if args.norm_method!='none':
      size_f=normalizeCounts(cttab_sel,method=args.norm_method,returnfactor=True,reversefactor=True,controlsgfile=args.control_sgrna)
    else:
      size_f=None
  else:
    size_f=normalizeCounts(cttab_sel,returnfactor=True,reversefactor=True)
  if size_f !=None:
    logging.info('size factor: '+','.join([str(x) for x in size_f]))
  
  # desmat=np.matrix([[1,1,1,1],[0,0,1,0],[0,0,0,1]]).getT()
  desmat=args.design_matrix
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    tginst.design_mat=desmat
  meanvardict={}
  for (tgid,tginst) in allgenedict.iteritems():
    #iteratenbem(tginst,debug=False,alpha_val=0.01,estimateeff=False,size_factor=size_f)
    ##sgrna_wide_dispersion_estimation_MAP_v2(tginst,tginst.design_mat)
    ngene+=1
    tginst.w_estimate=[]
    meanvardict[tgid]=tginst
    if ngene>maxfittinggene:
      break
  argsdict={'debug':False, 'alpha_val':0.01, 'estimateeff':False,'size_factor':size_f}
  runem_multiproc(meanvardict,args,nproc=args.threads,argsdict=argsdict)
  for (tgid,tginst) in meanvardict.iteritems():
    allgenedict[tgid]=tginst
  # model the mean and variance
  logging.info('Modeling the mean and variance ...')
  if maxfittinggene>0:
    mrm=MeanVarModel()
    # old: linear model
    mrm.get_mean_var_residule(allgenedict)
    mrm.model_mean_var_by_lm()
    # new method: generalized linear model
    #mrm.model_mean_disp_by_glm(allgenedict,args.output_prefix,size_f)
  else:
    mrm=None
  
  if returndict:
    return (allgenedict,mrm,size_f)
  # run the test again...
  logging.info('Run the algorithm for the second time ...')
  if hasattr(args,'threads') and args.threads>1:
    # multipel threads
    argsdict={'debug':False,'estimateeff':True,'meanvarmodel':mrm,'restart':False,'removeoutliers':args.remove_outliers,'size_factor':size_f,'updateeff':args.update_efficiency,'logem':False}
    runem_multiproc(allgenedict,args,nproc=args.threads,argsdict=argsdict)
  else:
    # only 1 thread
    ngene=0
    for (tgid,tginst) in allgenedict.iteritems():
      #try:
      if ngene % 1000 ==1 or args.debug:
        logging.info('Calculating '+tgid+' ('+str(ngene)+') ... ')
      if hasattr(args,'debug_gene') and args.debug_gene!=None and tginst.prefix != args.debug_gene:
        continue
      iteratenbem(tginst,debug=False,estimateeff=True,meanvarmodel=mrm,restart=False,removeoutliers=args.remove_outliers,size_factor=size_f,updateeff=args.update_efficiency)
      # Tracer()()
      ngene+=1
      if ngene>maxgene:
        break
      #except:
      #  logging.error('Error occurs while calculating beta values of gene '+tgid+'.')
      #  sys.exit(-1)
  # set up the w vector
  for (tgid,tginst) in allgenedict.iteritems():
    if len(tginst.w_estimate)==0:
      tginst.w_estimate=np.ones(len(tginst.sgrnaid))
  #Tracer()()
  # permutation
  iteratenbem_permutation(allgenedict,args,nround=args.permutation_round,removeoutliers=args.remove_outliers,size_factor=size_f)
  # correct for FDR
  from mageck.mleclassdef import gene_fdr_correction;
  gene_fdr_correction(allgenedict,args.adjust_method);
  # correct for CNV
  if args.cnv_norm is not None: # and args.cell_line is not None:
    logging.info('Performing copy number normalization.')
    (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(args.cnv_norm,args.beta_labels[1:])
    for i in range(len(args.beta_labels[1:])):
      if args.beta_labels[1:][i] not in CN_celldict:
        logging.warning(args.beta_labels[1:][i] + ' is not represented in the inputted copy number variation data.')
      else:
        logging.info('Normalizing by copy number with ' + args.beta_labels[1:][i] + ' as the reference cell line.')
    betascore_piecewisenorm(allgenedict,args.beta_labels,CN_arr,CN_celldict,CN_genedict)
  # write to file
  genefile=args.output_prefix+'.gene_summary.txt'
  sgrnafile=args.output_prefix+'.sgrna_summary.txt'
  logging.info('Writing gene results to '+genefile)
  logging.info('Writing sgRNA results to '+sgrnafile)
  write_gene_to_file(allgenedict,genefile,args,betalabels=args.beta_labels)
  write_sgrna_to_file(allgenedict,sgrnafile)
  return (allgenedict,mrm)


if __name__ == '__main__':
  try:
    mageckmle_main()
  except KeyboardInterrupt:
    sys.stderr.write("User interrupt me! ;-) Bye!\n")
    sys.exit(0)


