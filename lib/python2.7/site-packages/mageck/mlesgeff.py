'''
sgRNA efficiency related functions
'''


from __future__ import print_function
import sys
import numpy as np
from mageck.mleclassdef import *

import logging


def read_sgrna_eff(args):
  '''
  Read sgRNA efficiency score from file, and convert to initial prediction
  '''
  if args.sgrna_efficiency != None:
    # efficiency prediction
    nline=0
    sgrna_eff_dict={}
    sgscore_max=-1000000.0
    sgscore_min=10000000.0
    sgscore_minbound=-1
    sgscore_maxbound=0.3
    for line in open(args.sgrna_efficiency):
      nline+=1
      field=line.strip().split()
      if len(field)<= args.sgrna_eff_name_column or len(field)<=args.sgrna_eff_score_column:
        logging.warning('Not enough fields in line '+str(nline)+' of the sgRNA efficiency prediction file. Please check your --sgrna-eff-name-column and --sgrna-eff-score-column options.')
        continue
      sgid=field[args.sgrna_eff_name_column]
      try:
        sgscore=float(field[args.sgrna_eff_score_column])
      except ValueError:
        logging.warning('Error parsing sgRNA efficiency scores: '+field[args.sgrna_eff_score_column]+' at line '+str(nline)+' of the sgRNA efficiency prediction file. Skip this line ..')
        sgscore=None
      sgrna_eff_dict[sgid]=sgscore
      if sgscore > sgscore_max:
        sgscore_max=sgscore
      if sgscore < sgscore_min:
        sgscore_min=sgscore
    # end for
    logging.info(str(nline)+' lines processed in sgRNA efficiency prediction file '+args.sgrna_efficiency+'.')
    for (sgid,sgscore) in sgrna_eff_dict.iteritems():
      if sgscore == None:
        sgscore=0
      else:
        #sgscore= (sgscore-sgscore_min)*1.0/(sgscore_max-sgscore_min)
        sgscore = (sgscore-sgscore_minbound)*1.0/(sgscore_maxbound-sgscore_minbound)
      if sgscore < 1e-2:
        sgscore=1e-2
      if sgscore >1.0:
        sgscore=1.0
      sgrna_eff_dict[sgid]=sgscore
    args.sgrna_efficiency = sgrna_eff_dict


def sgrna_eff_initial_guess(args,allgenedict):
  '''
  Convert sgRNA efficiency to initial guess of P(eff)
  '''
  if args.sgrna_efficiency != None:
    logging.info('Converting sgRNA efficiency prediction to the initial guess of pi...')
    for (geneid,gk) in allgenedict.iteritems():
      sgid=gk.sgrnaid
      n=gk.nb_count.shape[1]
      gk.w_estimate=np.ones(n)
      for ii in range(len(sgid)):
        indsgid=sgid[ii]
        if indsgid in args.sgrna_efficiency:
          gk.w_estimate[ii]=args.sgrna_efficiency[indsgid]


