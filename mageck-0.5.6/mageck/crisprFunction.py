#!/usr/bin/env python
"""MAGeCK test module
Copyright (c) 2014 Wei Li, Han Xu, Xiaole Liu lab
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).
@status:  experimental
@version: $Revision$
@author:  Wei Li
@contact: li.david.wei AT gmail.com
"""


from __future__ import print_function
import sys
import math
import types
import logging

from mageck.mageckCount import *

from mageck.fileOps import *
from mageck.testVisual import *

from mageck.fdr_calculation import *




# debug
# try:
#   from IPython.core.debugger import Tracer
# except:
#   pass

def mmedian(lst):
  """
  get the median value
  """
  sortedLst = sorted(lst)
  lstLen = len(lst)
  if lstLen==0:
    return 0.0
  index = (lstLen - 1) // 2
  
  if (lstLen % 2):
    return sortedLst[index]
  else:
    return (sortedLst[index] + sortedLst[index + 1])/2.0

def getgeomean(v):
  meanval=sum([math.log(vx+0.1,2) for vx in v])/float(len(v))
  return 2**meanval-0.1

def getMeans(matt):
  # arithmatic mean
  #meanvalue=[sum(v)/float(len(v)) for v in matt]
  # geometric mean
  meanvalue=[getgeomean(v) for v in matt]
  return meanvalue

def getVars(matt):
  meanvalue=getMeans(matt)
  varvalue=[ sum([ (kj-meanvalue[i])*(kj-meanvalue[i])  for kj in matt[i]  ]  )/(float(len(matt[i]))-1)    for i in range(len(meanvalue))]
  #varvalue={k:sum([ (x-meanvalue[k])*(x-meanvalue[k]) for x in v])/(float(len(v))-1)  for (k,v) in ctable.iteritems()}
  return varvalue

def leastsquare(x,y,weight=None):
  """
  least squares fitting
  coefficients from y= a+bx
  return (b,a)
  reference: http://mathworld.wolfram.com/LeastSquaresFitting.html
  For weighted least square: http://goo.gl/pGpTZ6
  """
  n=len(x)
  if n != len(y):
    logging.error('Unequal length of vectors of x and y in least square')
    sys.exit(-1)
  if weight is None:
    sy=sum(y)
    sx=sum(x)
    sx2=sum([t*t for t in x])
    sxy=sum([x[i]*y[i] for i in range(n)])
    a=(sy*sx2-sx*sxy)/(n*sx2-sx*sx)
    b=(n*sxy-sx*sy)/(n*sx2-sx*sx)
    return (b,a)
  else:
    nw=sum(weight)
    sy=sum([y[i]*weight[i] for i in range(n)])
    sx=sum([x[i]*weight[i] for i in range(n)])
    sx2=sum([x[i]*x[i]*weight[i] for i in range(n)])
    sxy=sum([x[i]*y[i]*weight[i] for i in range(n)])
    a=(sy*sx2-sx*sxy)/(nw*sx2-sx*sx)
    b=(nw*sxy-sx*sy)/(nw*sx2-sx*sx)
    return (b,a)

def modelmeanvar(ctable,method='edger'):
  """
  model the relation between mean and variance
  """
  # calculate the mean and variance
  tablemat=ctable.values()
  meanvalue=getMeans(tablemat)
  varvalue=getVars(tablemat)
  # choose values with variance greater than mean
  meangood=[meanvalue[i] for i in range(len(meanvalue)) if meanvalue[i]<varvalue[i]]
  vargood=[varvalue[i]-meanvalue[i] for i in range(len(varvalue)) if meanvalue[i]<varvalue[i]]
  if len(meangood)<10:
    logging.warning('The variances between control replicates are too small. If they are technical replicates, merge these into one sample.')
    meangood=[meanvalue[i] for i in range(len(meanvalue)) ]
    vargood=[(lambda x: x if x>0.01 else 0.01 )(varvalue[i]-meanvalue[i]) for i in range(len(varvalue)) ]
  # log
  meanglog=[math.log(x+1,2) for x in meangood]
  varglog=[math.log(x+1,2) for x in vargood]
  # Tracer()()
  if method=='linear':
    # least square
    (k,b)=leastsquare(meanglog,varglog,meangood)
    if k<1:
      k=1
    if b<0:
      b=0
    return (k,b)
  elif method=='edger':
    dy=varglog
    dx=[2*x for x in meanglog]
    ret=(sum(dy)-sum(dx))*1.0/len(dx)
    return ret
  else:
    return 0


def getadjustvar(coef,meanval,method='mixed'):
  """
  From the model, get the adjusted variance
  """
  if method=='linear':
    k=coef[0];b=coef[1]
    if type(meanval) is types.FloatType:
      return (meanval**k)*(2**b)+meanval
    if type(meanval) is types.ListType:
      return [(z**k)*(2**b)+z for z in meanval]
  elif method=='edger':
    if type(meanval) is types.FloatType:
      return (meanval**2)*(2**coef)+meanval
    if type(meanval) is types.ListType:
      return [(z**2)*(2**coef)+z for z in meanval]
  elif method=='mixed':
    var1=getadjustvar(coef,meanval,method='linear')
    var2=getadjustvar(coef[2],meanval,method='edger')
    return [ (lambda x,y: x if x>y else y)(var1[i],var2[i]) for i in range(len(var1))]
  else:
    return meanval

def getnormcdf(x,lowertail=True):
  """
  Get the normal CDF function. used to calculate p-value
  """
  # ax=math.fabs(x)
  #axv=math.erfc(x/(2**0.5))/2; # higher tail
  if lowertail==False:
    #return axv
    return math.erfc(x/(2**0.5))/2
  else:
    #return 1-axv
    return math.erfc(-x/(2**0.5))/2
  #if (x>0 and lowertail==False) or (x<0 and lowertail==True):
  #  return axv
  #else:
  #  return 1-axv

def getNormalPValue(mean0,var0,mean1, lower=False):
  """
  Use truncated normal distribution to calculate the pvalue
  """
  # use ttmean to calculate the pvalue
  n=len(mean0)
  minmean1=min([x for x in mean1 if x>0])
  mean1_adj=[(lambda x: x if x >minmean1 else minmean1)(t) for t in mean1]
  # first, convert to standard normal distribution values
  t_theta=[(mean1_adj[i]-mean0[i])/math.sqrt(var0[i]) for i in range(n)]
  t_theta_0=[(0.0-mean0[i])/math.sqrt(var0[i]) for i in range(n)]
  #
  t_p=[getnormcdf(x,lowertail=lower) for x in t_theta]
  t_p_0=[getnormcdf(x,lowertail=True) for x in t_theta_0]
  if lower==True:
    return [(t_p[i]-t_p_0[i])/(1-t_p_0[i]) for i in range(n)]
  else:
    return [t_p[i]/(1-t_p_0[i]) for i in range(n)]


def getNBPValue(mean0,var0,mean1, lower=False,log=False):
  """
  Use negative binomial to calculate p-value
  Reference:
  http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html#scipy.stats.nbinom
  """
  from scipy.stats import nbinom
  n=len(mean0)
  nb_p=[mean0[i]/var0[i] for i in range(n)]; # consisitent with R
  nb_n0=[mean0[i]*mean0[i]/(var0[i]-mean0[i]) for i in range(n)]
  nb_n=[ (lambda t: t if t>=1 else 1)(x) for x in nb_n0]
  #
  if lower==True:
    if log==False:
      nb_p_low=nbinom.cdf(mean1,nb_n,nb_p)
    else:
      nb_p_low=nbinom.logcdf(mean1,nb_n,nb_p)
    return list(nb_p_low)
  else:
    if log==False:
      nb_p_low=nbinom.sf(mean1,nb_n,nb_p)
    else:
      nb_p_low=nbinom.logsf(mean1,nb_n,nb_p)
    return list(nb_p_low)

def calculate_gene_lfc(args,lfcval,sort_id,n_lower,sgrna2genelist,destkeys,ispos=False):
  """
  Calculate gene LFC using different methods
  Parameters:
    args
        Arguments
    lfcval
        sgRNA log fold change vector
    sortid
        sgRNA sort index
    n_lower
        alpha cutoff (integer)
    sgrna2genelist
        a {sgrnaid:gene} dict
    destkeys
        a [sgrnaid] vector
    ispos
        a boolean vector to indicate whether this is a positive selection
  Return value:
    genelfc
        a {geneid:lfc} dict
  """
  genesglfc={}
  ni=0
  for i in sort_id:
    ni+=1
    targetgene=sgrna2genelist[destkeys[i]]
    if targetgene not in genesglfc:
      genesglfc[targetgene]=[]
    if args.gene_lfc_method=='alphamean' or args.gene_lfc_method=='alphamedian':
      if ni*1.0<=n_lower:
        genesglfc[targetgene]+=[lfcval[i]]
    else:
      genesglfc[targetgene]+=[lfcval[i]]
  genelfc={}
  for (gid,vl) in genesglfc.items():
    if args.gene_lfc_method=='median' or args.gene_lfc_method=='alphamedian':
      lfc=mmedian(vl)
    elif args.gene_lfc_method=='secondbest':
      if ispos:
        vll=sorted(vl,reverse=True)
      else:
        vll=sorted(vl)
      if len(vll)>1:
        lfc=vll[1]
      else:
        lfc=0.0
    elif args.gene_lfc_method=='mean' or args.gene_lfc_method=='alphamean':
      if len(vl)>0:
        lfc=sum(vl)/len(vl)
      else:
        lfc=0.0
    else:
      lfc=0.0
    genelfc[gid]=lfc
  return genelfc

def crispr_test(tab,ctrlg,testg, destfile,sgrna2genelist,args):
  """
  main function of crispr test
  Parameters:
    tab
        Read count table
    ctrlg
        Index for control samples
    testg
        Index for treatment samples
    destfile
        Prefix for output file (sgrna_summary.txt)
    sgrna2genelist
        {sgrna:gene} mapping
    args
        Arguments
  Return value:
    (lowp,highp,sgrnalfc)
    lowp
        alpha cutoff for neg. selection
    highp
        alpha cutoff for pos. selection
    lower_gene_lfc
        {gene:lfc} dict. lfc is for neg. selection
    higher_gene_lfc
        {gene:lfc} dict. lfc is for pos. selection
  """
  n=len(tab)
  # control and test matrix
  tabctrl={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.iteritems()}
  tabtest={k:[v[i] for i in range(len(v)) if i in testg] for (k,v) in tab.iteritems()}
  # control matrix for mean-var estimation
  if len(ctrlg)>1 and args.variance_from_all_samples==False: # more than 1 controls
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in ctrlg] for (k,v) in tab.iteritems()}
  else: # only 1 control: use all the samples for estimation
    tabctrlmod={k:[v[i] for i in range(len(v)) if i in (ctrlg+testg)] for (k,v) in tab.iteritems()}
  # training using control samples
  model1=modelmeanvar(tabctrlmod,method='linear')
  #model2=modelmeanvar(tabctrl,method='edger')
  model=[x for x in model1];#+[model2]
  if type(model) is types.ListType:
    logging.debug('Adjusted model: '+'\t'.join([str(x) for x in model]))
  else:
    logging.debug('Adjusted model: k='+str(model))

  tabctrl_mat=tabctrl.values()
  tabctrlmodel_mat=tabctrlmod.values()
  tabc_mean=getMeans(tabctrl_mat)
  tabcmodel_mean=getMeans(tabctrlmodel_mat)
  #
  # setup the valid sgRNA flag
  validsgrna=[1]*n
  if hasattr(args,"remove_zero") and ( args.remove_zero=="control" or args.remove_zero=="both"):
    validsgrna=[ (lambda x: 1 if x>0 else 0)(t) for t in tabc_mean]
  # if mean of the control samples is 0: set it to greater than 0
  tabc_min=min([x for x in tabc_mean if x>0])
  tabc_mean=[ (lambda x: x if x>tabc_min else tabc_min)(t) for t in tabc_mean]
  tabc_var=getVars(tabctrlmodel_mat)
  tabc_adjvar=getadjustvar(model,tabc_mean,method='linear')

  # testing using tebtest
  nt=tabtest[tabtest.keys()[0]]
  ttmat=tabtest.values()
  ttmean=getMeans(ttmat)
  # set up the valid sgRNA flag
  if hasattr(args,"remove_zero") and ( args.remove_zero=="treatment" or args.remove_zero=="both"):
    validsgrna2=[ (lambda x: 1 if x>0 else 0)(t) for t in ttmean]
    validsgrna=[validsgrna[t]*validsgrna2[t] for t in range(n)]
  # use ttmean to calculate the pvalue
  # first, convert to standard normal distribution values
  tt_theta=[(ttmean[i]-tabc_mean[i])/math.sqrt(tabc_adjvar[i]) for i in range(n)]
  tt_abstheta=[math.fabs(tt_theta[i]) for i in range(n)]
  #
  try:
    # for consistency, use normal p values
    tt_p_lower=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=True)
    tt_p_higher=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=False)
    #tt_p_lower=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=True)
    #tt_p_higher=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=False)

    # tt_p_lower_score=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=True,log=True)
    # tt_p_higher_score=getNBPValue(tabc_mean,tabc_adjvar,ttmean,lower=False,log=True)
  #except ImportError:
  #  #logging.warning('An error occurs while trying to compute p values using scipy. Will use normal model instead of Negative Binomial model, but please check with your scipy installation.')
  #  #tt_p_lower=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=True)
  #  #tt_p_higher=getNormalPValue(tabc_mean,tabc_adjvar,ttmean,lower=False)
  except:
    logging.error('An error occurs while trying to compute p values. Quit..')
    sys.exit(-1)
  #
  #
  tt_p_twosided=[ (lambda x,y: 2*x if x<y else 2*y)(tt_p_lower[i],tt_p_higher[i]) for i in range(n)]
  tt_p_fdr=pFDR(tt_p_twosided,method=args.adjust_method)
  #
  # map sgRNA to genes
  gene_list = []
  sgrna_list = tabctrl.keys()
  for sgrna in sgrna_list:
    if sgrna2genelist is not None:
      gene_list.append(sgrna2genelist[sgrna])
    else:
      gene_list.append('NA')
  # normalize sgRNA scores and sort according to score
  CNVnorm = False
  if args.cnv_norm is not None and args.cell_line is not None:
    from mageck.cnv_normalization import read_CNVdata,sgRNAscore_piecewisenorm
    logging.info('Performing copy number normalization.')
    (CN_arr,CN_celldict,CN_genedict) = read_CNVdata(args.cnv_norm,[args.cell_line])
    if args.cell_line in CN_celldict:
      logging.info('Normalizing by copy number with' + args.cell_line + 'as the reference cell line.')
      CNVnorm = True
      norm_tt_theta = sgRNAscore_piecewisenorm(tt_theta,gene_list,CN_arr,CN_genedict)
      norm_tt_abstheta=[math.fabs(norm_tt_theta[i]) for i in range(n)]
      sort_id=[i[0] for i in sorted(enumerate(norm_tt_abstheta), key=lambda x:x[1],reverse=True)]
      # replace the original values of tt_theta
      tt_theta=norm_tt_theta
      tt_abstheta=norm_tt_abstheta
    else:
      logging.warning(args.cell_line + ' is not represented in the inputted copy number variation data.')
      sort_id=[i[0] for i in sorted(enumerate(tt_abstheta), key=lambda x:x[1],reverse=True)]
  else:
    sort_id=[i[0] for i in sorted(enumerate(tt_abstheta), key=lambda x:x[1],reverse=True)]
  #
  # lower_score and higher_score are used to sort sgRNAs
  tt_p_lower_score=tt_theta
  tt_p_higher_score=[-1*x for x in tt_theta]
  # write to file
  destfname=destfile+'.sgrna_summary.txt'
  destf=open(destfname,'w')
  destkeys=tabctrl.keys()
  dfmt="{:.5g}"
  
  # sgRNA log fold change 
  sgrnalfc=[0.0]*n
  # output to file
  header = ['sgrna','Gene','control_count','treatment_count','control_mean','treat_mean', 'LFC', 'control_var','adj_var','score','p.low','p.high','p.twosided','FDR','high_in_treatment']
  #if CNVnorm:
  #  header += ['CNVadj_score']
  print('\t'.join(header),file=destf)
  for i in sort_id:
    # sgRNA mapping to genes?
    if sgrna2genelist is not None:
      destkeygene=sgrna2genelist[destkeys[i]]
    else:
      destkeygene='None'
    report=[destkeys[i], destkeygene, '/'.join([dfmt.format(x) for x in tabctrl_mat[i]]), '/'.join([dfmt.format(x) for x in ttmat[i]])]
    t_r=[tabc_mean[i],ttmean[i]]
    lfcval=math.log(ttmean[i]+1.0,2)-math.log(tabc_mean[i]+1.0,2)
    t_r+=[ lfcval ] # log fold change
    sgrnalfc[i]=lfcval # save log fold change
    t_r+=[tabc_var[i],tabc_adjvar[i],tt_abstheta[i],tt_p_lower[i],tt_p_higher[i],tt_p_twosided[i],tt_p_fdr[i]]
    report+=[dfmt.format(x) for x in t_r]
    report+=[ttmean[i]>tabc_mean[i]]
    #if CNVnorm:
    #  report+=[dfmt.format(norm_tt_abstheta[i])] # add CNV-adjusted sgRNA scores
    print('\t'.join([str(x) for x in report]),file=destf)
  destf.close()
  #
  # prepare files for gene test
  if sgrna2genelist is not None:
    destfname=destfile+'.plow.txt'
    destkeys=tabctrl.keys()
    sort_id=[i[0] for i in sorted(enumerate(tt_p_lower_score), key=lambda x:x[1],reverse=False)]
    # output to file
    destf=open(destfname,'w')
    print('\t'.join(['sgrna','symbol','pool','p.low','prob','chosen']),file=destf)
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_lower_score[i], '1', validsgrna[i]]
      print('\t'.join([str(x) for x in report]),file=destf)
    destf.close()
    tt_p_lower_fdr=pFDR(tt_p_lower,method=args.adjust_method)
    n_lower=sum([1 for x in tt_p_lower if x <= args.gene_test_fdr_threshold])
    n_lower_p=n_lower*1.0/len(tt_p_lower)
    logging.debug('lower test FDR cutoff: '+str(n_lower_p))
    # calculate gene lfc
    lower_gene_lfc=calculate_gene_lfc(args,sgrnalfc,sort_id,n_lower,sgrna2genelist,destkeys)
    #
    destfname=destfile+'.phigh.txt'
    destf=open(destfname,'w')
    destkeys=tabctrl.keys()
    sort_id=[i[0] for i in sorted(enumerate(tt_p_higher_score), key=lambda x:x[1],reverse=False)]
    # output to file
    print('\t'.join(['sgrna','symbol','pool','p.high','prob','chosen']),file=destf)
    for i in sort_id:
      report=[destkeys[i], sgrna2genelist[destkeys[i]],'list', tt_p_higher_score[i], '1', validsgrna[i]]
      print('\t'.join([str(x) for x in report]),file=destf)
    destf.close()
    tt_p_higher_fdr=pFDR(tt_p_higher,method=args.adjust_method)
    n_higher=sum([1 for x in tt_p_higher if x <= args.gene_test_fdr_threshold])
    if n_higher>0:
      n_higher_p=n_higher*1.0/len(tt_p_higher)
    else:
      n_higher_p=0.01
    logging.debug('higher test FDR cutoff: '+str(n_higher_p))
    # calculate gene lfc
    higher_gene_lfc=calculate_gene_lfc(args,sgrnalfc,sort_id,n_higher,sgrna2genelist,destkeys,ispos=True)
    #
    return (n_lower_p,n_higher_p,lower_gene_lfc,higher_gene_lfc)
  else:
    return (None,None,None,None)

def rank_association_test(file,outfile,cutoff,args,adjustcutoff=True):
  if adjustcutoff: # adjust the alpha threshold to 0.05-0.5
    if cutoff<0.05:
      cutoff=0.05
    if cutoff>0.5:
      cutoff=0.5
  #rrapath='/'.join(sys.argv[0].split('/')[:-1]+["../bin/RRA"])
  rrapath='RRA'
  command=rrapath+" -i "+file+" -o "+outfile+" -p "+str(cutoff)
  if hasattr(args,'control_sgrna') and args.control_sgrna != None :
    command+=" --control "+args.control_sgrna
  if hasattr(args,'skip_gene'): 
    if args.skip_gene != None :
      for g in args.skip_gene:
        command+=" --skip-gene "+g
    else:
      command+=" --skip-gene NA --skip-gene na "
  else:
    command+=" --skip-gene NA "
  # command+=" --min-number-goodsgrna 2 "
  if hasattr(args,"additional_rra_parameters") and args.additional_rra_parameters != None:
    command+=" "+args.additional_rra_parameters+" "
  systemcall(command)


def magecktest_removetmp(prefix):
  tmpfile=[prefix+'.plow.txt',prefix+'.phigh.txt',prefix+'.gene.low.txt',prefix+'.gene.high.txt']
  for f in tmpfile:
    systemcall('rm '+f,cmsg=False)


def magecktest_parsetreatmentsfromday0(args,samplelabelindex):
  """
  Reconstruct the groups of treatment and control from --day0-label
  """
  samples=[s for s in samplelabelindex.keys()]
  day0labelstr=args.day0_label
  args.day0_label=args.day0_label.split(',')
  for dl in args.day0_label:
    if dl not in samples:
      logging.error('Label '+dl+' specified in the --day0-label option does not match count table. Please double check.')
      sys.exit(-1)
  nonday0sample=[x for x in samples if x not in args.day0_label]
  if len(nonday0sample)==0:
    logging.error('At least 1 non day0-label sample should be specified.')
    sys.exit(-1)
  args.treatment_id=nonday0sample
  args.control_id=[day0labelstr]*len(nonday0sample)


def magecktest_main(args):
  """
  Main entry for MAGeCK test function
  """

  # stat test
  if args.subcmd == 'run' or args.subcmd == 'test':
    # read counts from file
    if args.subcmd == 'test':
      mapres=getcounttablefromfile(args.count_table)
    else:
      mapres=getcounttablefromfile(args.output_prefix+'.count.txt')
    cttab=mapres[0]
    sgrna2genelist=mapres[1]
    samplelabelindex=mapres[2]

    if len(cttab)==0:
      sys.exit(-1)
    nsample=len(cttab[cttab.keys()[0]])

    # process day0-label
    if args.day0_label != None:
      magecktest_parsetreatmentsfromday0(args,samplelabelindex)

    # iterate control group and treatment group
    supergroup_control=args.control_id
    supergroup_treat=args.treatment_id
    # control group and treatment group labels
    labellist_control=[]
    labellist_treat=[]
    # R visualization init
    vrv=VisualRValue()
    vrv.outprefix=args.output_prefix
    vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt'
    vrv.startRTemplate()
    vrvrnwcplabel=[]; # labels to write in rnw

    # loop by comparisons
    for cpindex in range(len(supergroup_treat)):
      # convert the sample label to sample index
      if cpindex==0:
        cp_prefix=args.output_prefix
      else:
        cp_prefix=args.output_prefix+'.'+str(cpindex)
      # labels
      (treatgroup,treatgrouplabellist)=parse_sampleids(supergroup_treat[cpindex],samplelabelindex)
      treatgroup_label=str(supergroup_treat[cpindex])
      logging.info('Treatment samples:'+treatgroup_label)
      logging.info('Treatment sample index:'+','.join([str(x) for x in treatgroup]))
      labellist_treat+=[treatgroup_label]
      if supergroup_control != None:
        (controlgroup,controlgrouplabellist)=parse_sampleids(supergroup_control[cpindex],samplelabelindex)
        controlgroup_label=str(supergroup_control[cpindex]); # only for display
        logging.info('Control samples:'+controlgroup_label)
      else:
        #controlgroup=[x for x in range(nsample) if x not in treatgroup]
        #controlgrouplabellist=[samplelabelindex[x] for x in range(nsample) if x not in treatgroup]
        xls=[x for x in range(nsample) if x not in treatgroup]
        (controlgroup,controlgrouplabellist)=parse_sampleids(','.join([str(t) for t in xls]),samplelabelindex)
        controlgroup_label='rest'
        logging.info('Control samples: the rest of the samples')
      logging.info('Control sample index:'+','.join([str(x) for x in controlgroup]))
      labellist_control+=[controlgroup_label]
      # read the sgRNA-gene table for rank association
      # normalization
      cttab_sel={k:([v[i] for i in controlgroup + treatgroup]) for (k,v) in cttab.iteritems()}; # controlgroup do not overlap with treatgroup
      if hasattr(args,'norm_method'):
        nttab=normalizeCounts(cttab_sel,method=args.norm_method,controlsgfile=args.control_sgrna)
      else:
        nttab=normalizeCounts(cttab_sel)
      # write normalized counts to file
      if hasattr(args,'normcounts_to_file'):
        if args.normcounts_to_file:
          # counts
          mageck_printdict(nttab,args,sgrna2genelist,samplelabelindex,controlgroup+treatgroup)
      controlgroup_ids=list(range(len(controlgroup)))
      treatgroup_ids=list(range(len(controlgroup),len(controlgroup+treatgroup)))
      # perform sgRNA test, and prepare files for gene test
      gene_as_cutoff=crispr_test(nttab, controlgroup_ids, treatgroup_ids, cp_prefix,sgrna2genelist,args)
      #
      if gene_as_cutoff[0] is not None:
        rank_association_test(cp_prefix+'.plow.txt',cp_prefix+'.gene.low.txt',gene_as_cutoff[0],args)
      if gene_as_cutoff[1] is not None:
        rank_association_test(cp_prefix+'.phigh.txt',cp_prefix+'.gene.high.txt',gene_as_cutoff[1],args,adjustcutoff=False) # update: fpr positive selection, do not adjust alpha cutoff
      # merge different files
      merge_rank_files(cp_prefix+'.gene.low.txt',cp_prefix+'.gene.high.txt',cp_prefix+'.gene_summary.txt',args,gene_as_cutoff)
      if cpindex>0:
        if cpindex>1:
          label1=''
        else:
          if len(labellist_treat)>0:
            label1=labellist_treat[0]+'_vs_'+labellist_control[0]+'|'
          else:
            label1=''
        label2=treatgroup_label+'_vs_'+controlgroup_label+'|'
        merge_rank_summary_files(args.output_prefix+'.gene_summary.txt',cp_prefix+'.gene_summary.txt',args.output_prefix+'.gene_summary.txt',args,lowfile_prefix=label1,highfile_prefix=label2)
      # visualization: load top k genes
      # print(str(samplelabelindex))
      vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' neg.'
      vrvrnwcplabel+=[vrv.cplabel]
      vrv.cpindex=[2+12*cpindex+1]
      vrv.loadTopKWithExp(cp_prefix+'.gene.low.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist)
      vrv.cplabel=treatgroup_label+'_vs_'+controlgroup_label+' pos.'
      vrvrnwcplabel+=[vrv.cplabel]
      vrv.cpindex=[2+12*cpindex+6+1]
      vrv.loadTopKWithExp(cp_prefix+'.gene.high.txt',nttab,sgrna2genelist,controlgrouplabellist+treatgrouplabellist)

      # clean the file
      if args.keep_tmp==False:
        magecktest_removetmp(cp_prefix)
        if cpindex>0:
          systemcall('rm '+cp_prefix+'.gene_summary.txt',cmsg=False)
          systemcall('rm '+cp_prefix+'.sgrna_summary.txt',cmsg=False)
      # end cleaning
    # end cpindex loop

    # generate pdf file
    # write to rnw file buffer
    vrv.genesummaryfile=args.output_prefix+'.gene_summary.txt'
    vrv.getGeneSummaryStat(args,isplot=False)
    vrv.comparisonlabel=vrvrnwcplabel; # replace the label field
    vrv.writeGeneSummaryStatToBuffer()
    # write to rnw and R file
    vrv.closeRTemplate()
    if hasattr(args, "pdf_report") and args.pdf_report:
      vrv.generatePDF(args.keep_tmp)
  # end if
