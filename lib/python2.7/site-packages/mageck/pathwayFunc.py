#!/usr/bin/env python
"""MAGeCK pathway analysis module
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
from mageck.crisprFunction import *
import random
import sys
import logging

from mageck.fileOps import *


def mageck_read_GMT(args):
  '''
  Read the GMT pathway files
  '''
  logging.info('Reading gmt:'+args.gmt_file)
  pathwaydict={}
  for line in open(args.gmt_file):
    field=line.strip().split()
    if len(field)<3:
      continue
    pathwaydict[field[0]]=[x.upper() for x in field[2:]]
  logging.info(str(len(pathwaydict))+' pathways loaded.')
  return pathwaydict

def mageck_readgeneranking(fname,args,log2=False,columnid=2):
  """
  Read gene ranking files into a dictionary

  Arguments
  ---------
  fname: string
    gene ranking file.
  args: argparse object
  log2:
    whether to log-2 transform the score
  columnid: int/string
    The column id of the score, default is 2 (3rd score)
  """
  geneinfo={}
  nline=0
  columnid_int=-1
  try:
    columnid_int=int(columnid)
  except ValueError:
    pass
  for line in open(fname):
    nline=nline+1
    field=line.strip().split()
    if nline<2: # skip the first line
      for i in range(len(field)):
        if field[i]==columnid:
          columnid_int=i
          break
      if columnid_int==-1:
        logging.error('Cannot determine the column to be used for ranking. Please double check with the id you provide.')
        sys.exit(-1)
      else:
        logging.info('Column used for ranking: '+str(columnid)+' ('+str(columnid_int)+')')
      continue
    if len(field)<columnid_int+1:
      logging.error('Cannot read field '+str(columnid_int)+' in file '+fname+', line '+str(nline))
      sys.exit(-1)
    genename=field[0].upper()
    gscore=0
    try:
      nf=float(field[columnid_int])
      if log2:
        gscore=math.log(nf,2)
      else:
        gscore=nf
    except ValueError:
      gscore=1
    geneinfo[genename]=gscore
  return geneinfo

def mageck_removetmprra(args):
  if args.single_ranking:
    tmpfile=[args.output_prefix+'.pathway.tmp']
  else:
    tmpfile=[args.output_prefix+'.pathway.low.tmp',args.output_prefix+'.pathway.high.tmp',args.output_prefix+'.pathway.low.txt',args.output_prefix+'.pathway.high.txt']

  for f in tmpfile:
    systemcall('rm '+f,cmsg=False)

def mageck_pathwayrra(args):
  """perform pathway anaylsis using RRA
  """
  # reading pathway files
  pdict=mageck_read_GMT(args)
  fname=args.gene_ranking
  if args.single_ranking:
    columnid=args.ranking_column
    mageck_pathwayrra_onedir(args,pdict,columnid,fname,args.output_prefix+'.pathway.tmp',args.output_prefix+'.pathway_summary.txt')
    # columnid=2; # default: 3rd column (neg. selected p values)
  else:
    tmppath_low=args.output_prefix+'.pathway.low.tmp'
    tmppath_high=args.output_prefix+'.pathway.high.tmp'
    rraout_low=args.output_prefix+'.pathway.low.txt'
    rraout_high=args.output_prefix+'.pathway.high.txt'
    columnid=args.ranking_column
    mageck_pathwayrra_onedir(args,pdict,columnid,fname,tmppath_low,rraout_low)
    columnid=args.ranking_column_2; # columnid=6 if sgRNA number in positive selection is not omitted
    mageck_pathwayrra_onedir(args,pdict,columnid,fname,tmppath_high,rraout_high)
    # merge different files
    cutoffinfo=(args.pathway_alpha,args.pathway_alpha,{},{})
    merge_rank_files(rraout_low,rraout_high,args.output_prefix+'.pathway_summary.txt',args,cutoffinfo)

  if args.keep_tmp==False:
    mageck_removetmprra(args)

def mageck_pathwayrra_onedir(args,pdict,cid,sourcefile,rra_path_input_file,rra_path_output_file):
  '''
  Calling RRA for pathway test
  '''
  logging.debug('Performing pathway ranking test; input:'+rra_path_input_file+', output:'+rra_path_output_file)
  # preparing input files
  # rra_path_input_file=args.output_prefix+'.pathway.tmp'
  rra_path_input=open(rra_path_input_file,'w')
  print('\t'.join(['gene','pathway','pool','score']),file=rra_path_input)
  
  ginfo=mageck_readgeneranking(sourcefile,args,columnid=cid)
  ginfo_st=sorted(ginfo.items(),key=lambda x: x[1])
  nthreshold=0
  for gtuple in ginfo_st:
    genename=gtuple[0]
    gscore=0
    try:
      nf=gtuple[1]
      if nf<0.05:
        nthreshold+=1
      gscore=math.log(nf,2)
    except ValueError:
      gscore=0
    ## old method: for each gene, make a copy for each pathway
    #niters=0
    #for (k,v) in pdict.items():
    #  if genename in v:
    #    niters+=1
    #    print('\t'.join([genename+'_copy'+str(niters),k,'list',str(gscore)]),file=rra_path_input)
    #if niters==0:
    #  print('\t'.join([genename,'NA','list',str(gscore)]),file=rra_path_input)
    ## new method: pathway name separated by ','
    pnamelist=[]
    for (k,v) in pdict.items():
      if genename in v:
        pnamelist+=[k]
    if len(pnamelist)==0:
      print('\t'.join([genename,'NA','list',str(gscore)]),file=rra_path_input)
    else:
      print('\t'.join([genename,','.join(pnamelist),'list',str(gscore)]),file=rra_path_input)
  rra_path_input.close()
  # rra_threshold=nthreshold*1.0/len(ginfo_st)
  rra_threshold=args.pathway_alpha

  # rank association test
  rank_association_test(rra_path_input_file,rra_path_output_file,rra_threshold,args)

def mageck_pathway_standardize(gdict):
  '''
  Standardize the scores in a dictionary
  '''
  ginfo_vec=sorted(gdict.values())
  ginfo_med=ginfo_vec[len(ginfo_vec)//2]
  ginfo_var=sum([ (x-ginfo_med)**2 for x in ginfo_vec])/(len(ginfo_vec)-1)
  ginfo_std=ginfo_var**0.5
  ginfo_rt={ k:((v-ginfo_med)/ginfo_std) for (k,v) in gdict.items()}
  return ginfo_rt

def mageck_pathway_ztest_permutation(args,gdict,pdict,pdictpvals):
  '''
  Perform permutation test
  '''
  pdictqvals={}
  gdictvals=gdict.values()
  niter=10000
  print('Permutation test, please wait',file=sys.stderr,end='')
  itercounter=0
  for iter in range(niter):
    itercounter+=1
    if itercounter % 100 ==1:
      print('.',end='',file=sys.stderr)
    for (pname,pitem) in pdictpvals.items():
      nsize=pitem[2]
      # testvecsq=[x**2 for x in testvec]
      if nsize>1:
        #testval2=sum(random.sample(gdictvals,nsize))/(nsize**0.5); # one direction scale
        testval2=sum([x**2 for x in random.sample(gdictvals,nsize)])/2/(nsize-1)-0.5; # two directional scale
      else:
        testval2=0
      if pname not in pdictqvals:
        pdictqvals[pname]=[0.0,0.0]
      if pitem[3]>=testval2:
        pdictqvals[pname][0]+=1
      if pitem[3]<=testval2:
        pdictqvals[pname][1]+=1
  print('',file=sys.stderr)
  pdictqvals={k:(v[0]*1.0/niter,v[1]*1.0/niter) for (k,v) in pdictqvals.items()}
  return pdictqvals




def mageck_pathway_ztest(args,gdict,pdict):
  '''
  Perform z test on the normalized gene score (gdict) and pathway (pdict)
  Return values:
    pdict pvals: a dictionary structure, with key the pathway name and value a n-tuple including:
                   p value (low), pvalue (high), size of the pathway genes in the destination, and test statistics
  '''
  pdictpvals={}
  for (pname,pitem) in pdict.items():
    testvec=[gdict[x] for x in pitem if x in gdict]
    testvecsq=[x**2 for x in testvec]
    if len(testvec)>1:
      #testval2=sum(testvec)/(len(testvec)**0.5); # one direction scale
      testval2=sum(testvecsq)/2/(len(testvec)-1)-0.5; # two directional scale
    else:
      testval2=0
    sq_p_low=getnormcdf(testval2)
    sq_p_high=getnormcdf(testval2,lowertail=False)
    pdictpvals[pname]=(sq_p_low,sq_p_high,len(testvec),testval2)
  return pdictpvals


def mageck_pathwaygsa(args):
  '''
  pathway enrichment analysis
  '''
  # reading pathway files
  pdict=mageck_read_GMT(args)
  # reading gene ranking files
  ginfo=mageck_readgeneranking(args.gene_ranking,args,log2=True)
  if args.gene_ranking_2 is not None:
    ginfo2=mageck_readgeneranking(args.gene_ranking_2,args,log2=True)
    # merge the scores of ginfo2 into ginfo
    for gk in ginfo.keys():
      if gk in ginfo2.keys():
        ginfo[gk]=ginfo[gk]-ginfo2[gk]
  # rank transform
  ranktransform=True
  ginfo_rt=ginfo
  if ranktransform:
    ginfo_st=sorted(ginfo.items(),key=lambda x: x[1])
    gauss_x=sorted([ random.gauss(0,1) for x in ginfo_st])
    ginfo_rt={ ginfo_st[i][0]:gauss_x[i] for i in range(len(ginfo_st))}
  # standardize the scores
  ginfo_sd=mageck_pathway_standardize(ginfo_rt)

  # test for pathways
  pathway_pval=mageck_pathway_ztest(args,ginfo_sd,pdict)
  #pathway_pval_tup=sorted(pathway_pval.items(),key=lambda x : min(x[1][:2]))
  #x_pvalues=[min(t[1][:2]) for t in pathway_pval_tup]
  #x_fdr=pFDR(x_pvalues)

  # permutation?
  pathway_qval=mageck_pathway_ztest_permutation(args,ginfo_sd,pdict,pathway_pval)
  pathway_qval_tup=sorted(pathway_qval.items(),key=lambda x:min(x[1]))
  x_qvalues=[min(t[1]) for t in pathway_qval_tup]
  x_fdr=pFDR(x_qvalues)
  # output file
  rsa_path_output_file=args.output_prefix+'.pathway_summary.txt'
  rsa_file=open(rsa_path_output_file,'w')
  # print('\t'.join(['PATHWAY','low','high','size','score','pvalue','q.low','q.high','FDR']),file=rsa_file)
  print('\t'.join(['PATHWAY','size','score','q.low','q.high','FDR']),file=rsa_file)
  for i in range(len(pathway_qval_tup)):
    lx=pathway_qval_tup[i]
    lpadj=x_fdr[i]
    print('\t'.join([lx[0],'\t'.join([str(z) for z in pathway_pval[lx[0]][2:]]), # str(min(lx[1][:2])),
        str(pathway_qval[lx[0]][0]), str(pathway_qval[lx[0]][1]), str(lpadj)]),file=rsa_file)
  rsa_file.close()

def mageck_pathwaygsa_fast(args):
  """
  Perform fast GSEA analysis using mageckGSEA
  """
  if args.single_ranking:
    columnid=args.ranking_column
    outputfile=args.output_prefix+'.pathway_summary.txt'
    gseacommand="mageckGSEA "
    gseacommand+=" -s "
    gseacommand+=" -c "+str(columnid)+" "
    gseacommand+=" -p "+str(args.permutation)
    gseacommand+=" -g \""+args.gmt_file+"\" "
    gseacommand+=" -r \""+args.gene_ranking+"\" "
    gseacommand+=" -o "+outputfile+" "
    logging.info("GSEA command: "+gseacommand)
    systemcall(gseacommand)
    # columnid=2; # default: 3rd column (neg. selected p values)
  else:
    rraout_low=args.output_prefix+'.pathway.low.txt'
    rraout_high=args.output_prefix+'.pathway.high.txt'
    #
    columnid=args.ranking_column
    gseacommand="mageckGSEA "
    gseacommand+=" -s "
    gseacommand+=" -c "+str(columnid)+" "
    gseacommand+=" -p "+str(args.permutation)
    gseacommand+=" -g \""+args.gmt_file+"\" "
    gseacommand+=" -r \""+args.gene_ranking+"\" "
    gseacommand+=" -o "+rraout_low+" "
    logging.info("GSEA command for negative selection: "+gseacommand)
    systemcall(gseacommand)
    #
    columnid=args.ranking_column_2; # columnid=6 if sgRNA number in positive selection is not omitted
    gseacommand="mageckGSEA "
    gseacommand+=" -s "
    gseacommand+=" -c "+str(columnid)+" "
    gseacommand+=" -p "+str(args.permutation)
    gseacommand+=" -g \""+args.gmt_file+"\" "
    gseacommand+=" -r \""+args.gene_ranking+"\" "
    gseacommand+=" -o "+rraout_high+" "
    logging.info("GSEA command for positive selection: "+gseacommand)
    systemcall(gseacommand)
    
    # merge different files
    merge_gsea_rank_files(rraout_low,rraout_high,args.output_prefix+'.pathway_summary.txt',args)
  if args.keep_tmp==False:
    mageck_removetmprra(args)

def merge_gsea_rank_files(lowfile,highfile,outfile,args):
  """
  Merge ranked file in GSEA
  """
  # process low_file
  nline=0
  gdict_low={}
  gdict_low_order=[]
  for line in open(lowfile):
    field=line.strip().split()
    nline+=1
    if nline==1: # skip the first line
      continue
    if len(field)<4:
      logging.error('The number of fields in file '+lowfile+' is <4.')
      sys.exit(-1)
    gid=field[0]
    gdict_low[gid]=field
    gdict_low_order+=[gid]
  # process high_file
  maxnline=nline
  gdict_high={}
  nline=0
  for line in open(highfile):
    field=line.strip().split()
    nline+=1
    if nline==1: # skip the first line
      continue
    if len(field)<4:
      logging.error('The number of fields in file '+lowfile+' is <4.')
      sys.exit(-1)
    gid=field[0]
    gdict_high[gid]=field
    if gid not in gdict_low:
      logging.warning('Item '+gid+' appears in '+highfile+', but not in '+lowfile+'.')
      gdict_low[gid]=[[field[0],field[1],0.0,1.0,1.0,1.0,maxnline,0,0.0]] # 
  # output
  ofhd=open(outfile,'w')
  print('\t'.join(['id','num','neg|score', 'neg|rra', 'neg|p-value','neg|fdr','neg|rank','neg|goodgene', 'neg|lfc', 'pos|score', 'pos|rra', 'pos|p-value','pos|fdr','pos|rank','pos|goodgene','pos|lfc']),file=ofhd)
  for gid in gdict_low_order:
    lowvecp=gdict_low[gid]
    if gid in gdict_high:
      lowvecp+=gdict_high[gid][2:]
    else:
      logging.warning('Item '+gid+' appears in '+lowfile+', but not in '+highfile+'.')
      lowvecp+=[0.0,1.0,1.0,1.0,maxnline,0,0.0] # 
    print('\t'.join(lowvecp),file=ofhd)
  ofhd.close()



def mageck_pathwaytest(args):
  '''
  The main entry for pathway test
  '''
  if args.method == 'gsea':
    mageck_pathwaygsa_fast(args)
  elif args.method == 'rra':
    mageck_pathwayrra(args)
  #mageck_pathwaygsa(args)
