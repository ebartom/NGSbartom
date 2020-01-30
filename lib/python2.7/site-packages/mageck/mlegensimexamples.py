#!/usr/bin/env python

from __future__ import print_function

import re
import sys
import scipy
from scipy.stats import nbinom
import random
import math
import numpy as np
import numpy.linalg as linalg

# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

# importing predefined functions
from mageck.mleclassdef import *
from mageck.mledesignmat import *
from mageck.mleem import *
from mageck.mleinstanceio import *
from mageck.mlemeanvar import *

from mageck.mlemageck import *



##
# simulated test case

def gentestcase1(nsg=10):
  '''
  The first testcase, 2 samples, control and treatment
  '''
  vark=0.01
  # basic parameters
  sks=SimCaseSimple()
  sks.beta0=[random.uniform(3,10) for i in range(nsg)] # these are the base 
  sks.beta1=[random.random()*5] # treatment
  print('beta_0:'+'\t'.join([decformat(x) for x in sks.beta0]))
  print('beta_1:'+decformat(sks.beta1[0]))
  # mean and variance 
  sks.mu=[[math.exp(t) for t in sks.beta0]]
  for t in sks.beta0:
    sks.mu+=[[math.exp(t+sks.beta1[0]) for t in sks.beta0]]
  #sks.var0=[t+vark*(t*t) for t in sks.mu0]
  #sks.var1=[[t+vark*(t*t) for t in sks.mu1[0]]]
  sks.var=[[t+vark*(t*t) for t in sks.mu[i]] for i in range(2)]
  #print('mu_0:'+'\t'.join([decformat(x) for x in sks.mu0]))
  #print('var_0:'+'\t'.join([decformat(x) for x in sks.var0]))
  #print('mu_1:'+'\t'.join([decformat(x) for x in sks.mu1[0]]))
  #print('var_1:'+'\t'.join([decformat(x) for x in sks.var1[0]]))
  # parameters for generating NB counts
  #sks.nb_p0=[sks.mu0[i]/sks.var0[i] for i in range(nsg)]
  #sks.nb_p1=[[sks.mu1[0][i]/sks.var1[0][i] for i in range(nsg)]]
  sks.nb_p=[[sks.mu[j][i]/sks.var[j][i] for i in range(nsg)] for j in range(2)]
  #sks.nb_r0=[sks.mu0[i]*sks.mu0[i]/(sks.var0[i]-sks.mu0[i]) for i in range(nsg)]
  #sks.nb_r1=[[sks.mu1[0][i]*sks.mu1[0][i]/(sks.var1[0][i]-sks.mu1[0][i]) for i in range(nsg)]]
  sks.nb_r=[[sks.mu[j][i]*sks.mu[j][i]/(sks.var[j][i]-sks.mu[j][i]) for i in range(nsg)] for j in range(2)]
  # 
  #sks.nb_count0=[nbinom.rvs(sks.nb_r0[i],sks.nb_p0[i]) for i in range(nsg)]
  #sks.nb_count1=[[nbinom.rvs(sks.nb_r1[0][i],sks.nb_p1[0][i]) for i in range(nsg)]]
  sks.nb_count=[[nbinom.rvs(sks.nb_r[j][i],sks.nb_p[j][i]) for i in range(nsg)] for j in range(2)]
  # design matrix
  # sks.design_mat=getsimpledesignmat(nsg)
  sks.design_mat=np.matrix([[1]])
  
  return (sks)

def gentestcase2(nsg=10):
  '''
  The second testcase, 2 samples, control and treatment
  '''
  vark=0.01
  # desmat=np.matrix([[0,0],[0,1],[1,0],[1,1]])
  desmat=np.matrix([[1,0,0],[0,1,0],[0,1,1],[1,1,1]])
  (nsample,nbeta)=desmat.shape
  # basic parameters
  sks=SimCaseSimple()
  sks.prefix='sample2'
  sks.design_mat=desmat
  sks.beta0=[random.uniform(3,10) for i in range(nsg)] # these are the base 
  sks.beta1=[random.random()*5 for i in range(nbeta)] # treatments;size: nbeta
  print('beta_0:'+'\t'.join([decformat(x) for x in sks.beta0]))
  print('beta_1:'+'\t'.join([decformat(x) for x in sks.beta1]))
  # mean and variance 
  mu0=[math.exp(t) for t in sks.beta0] # size: nsg
  tprod=desmat*np.matrix(sks.beta1).getT() # size: nsample*1  
  tprodlist=[x[0] for x in tprod.tolist()] # size: nsample*1  
  sks.mu=[mu0]
  for nr in range(nsample):
    sgi=[math.exp(t+tprodlist[nr]) for t in sks.beta0]
    sks.mu+=[sgi]
  # sks.var0=[t+vark*(t*t) for t in sks.mu0]
  sks.var=[[t+vark*(t*t) for t in tl] for tl in sks.mu]
  for i in range(nsample+1): # including 1 base and n samples
    print('mu_:'+str(i)+'\t'.join([decformat(x) for x in sks.mu[i]]))
    print('var_:'+str(i)+'\t'.join([decformat(x) for x in sks.var[i]]))
  # parameters for generating NB counts
  #sks.nb_p0=[sks.mu0[i]/sks.var0[i] for i in range(nsg)]
  #sks.nb_r0=[sks.mu0[i]*sks.mu0[i]/(sks.var0[i]-sks.mu0[i]) for i in range(nsg)]
  #sks.nb_p1=[[sks.mu1[t][i]/sks.var1[t][i] for i in range(nsg)] for t in range(nsample)]
  #sks.nb_r1=[[sks.mu1[t][i]*sks.mu1[t][i]/(sks.var1[t][i]-sks.mu1[t][i]) for i in range(nsg)] for t in range(nsample)]
  sks.nb_p=[[sks.mu[t][i]/sks.var[t][i] for i in range(nsg)] for t in range(nsample+1)]
  sks.nb_r=[[sks.mu[t][i]*sks.mu[t][i]/(sks.var[t][i]-sks.mu[t][i]) for i in range(nsg)] for t in range(nsample+1)]
  # 
  #sks.nb_count0=[nbinom.rvs(sks.nb_r0[i],sks.nb_p0[i]) for i in range(nsg)]
  #sks.nb_count1=[[nbinom.rvs(sks.nb_r1[t][i],sks.nb_p1[t][i]) for i in range(nsg)] for t in range(nsample)]
  sks.nb_count=[[nbinom.rvs(sks.nb_r[t][i],sks.nb_p[t][i]) for i in range(nsg)] for t in range(nsample+1)]
  
  return (sks)

def gentestcase3(nsg=10,desmat=None):
  '''
  The third testcase, with efficient 
  '''
  vark=0.01
  effiprob=0.5 # the probability that a sgRNA is efficient
  # desmat=np.matrix([[0,0],[0,1],[1,0],[1,1]])
  if desmat==None:
    # desmat=np.matrix([[1,0,0],[0,1,0],[0,1,1],[1,1,1]])
    desmat=np.matrix([[1,0,0,1],[0,1,1,1],[1,0,1,0]]).getT()
  (nsample,nbeta)=desmat.shape
  # basic parameters
  sks=SimCaseSimple()
  sks.prefix='sample3'
  sks.design_mat=desmat
  #sks.beta0=[random.uniform(3,10) for i in range(nsg)] # these are the base 
  #sks.beta1=[(random.random())*5 for i in range(nbeta)] # treatments;size: nbeta
  sks.beta0=[random.uniform(5,10) for i in range(nsg)] # these are the base 
  sks.beta1=[(random.random()*2-1)*5 for i in range(nbeta)] # treatments;size: nbeta
  print('beta_0:'+'\t'.join([decformat(x) for x in sks.beta0]))
  print('beta_1:'+'\t'.join([decformat(x) for x in sks.beta1]))
  # efficiency
  sks.isefficient=[ (lambda x: 1 if x>=effiprob else 0)(random.random()) for i in range(nsg)]
  # mean and variance 
  mu0=[math.exp(t) for t in sks.beta0] # size: nsg
  tprod=desmat*np.matrix(sks.beta1).getT() # size: nsample*1  
  tprodlist=[x[0] for x in tprod.tolist()] # size: nsample*1  
  sks.mu=[mu0]
  for nr in range(nsample):
    sgi=[math.exp(sks.beta0[ti]+tprodlist[nr]*sks.isefficient[ti]) for ti in range(nsg)]
    sks.mu+=[sgi]
  # sks.var0=[t+vark*(t*t) for t in sks.mu0]
  sks.var=[[t+vark*(t*t) for t in tl] for tl in sks.mu]
  for i in range(nsample+1): # including 1 base and n samples
    print('mu_:'+str(i)+'\t'.join([decformat(x) for x in sks.mu[i]]))
    print('var_:'+str(i)+'\t'.join([decformat(x) for x in sks.var[i]]))
  # parameters for generating NB counts
  #sks.nb_p0=[sks.mu0[i]/sks.var0[i] for i in range(nsg)]
  #sks.nb_r0=[sks.mu0[i]*sks.mu0[i]/(sks.var0[i]-sks.mu0[i]) for i in range(nsg)]
  #sks.nb_p1=[[sks.mu1[t][i]/sks.var1[t][i] for i in range(nsg)] for t in range(nsample)]
  #sks.nb_r1=[[sks.mu1[t][i]*sks.mu1[t][i]/(sks.var1[t][i]-sks.mu1[t][i]) for i in range(nsg)] for t in range(nsample)]
  sks.nb_p=[[sks.mu[t][i]/sks.var[t][i] for i in range(nsg)] for t in range(nsample+1)]
  sks.nb_r=[[sks.mu[t][i]*sks.mu[t][i]/(sks.var[t][i]-sks.mu[t][i]) for i in range(nsg)] for t in range(nsample+1)]
  # 
  #sks.nb_count0=[nbinom.rvs(sks.nb_r0[i],sks.nb_p0[i]) for i in range(nsg)]
  #sks.nb_count1=[[nbinom.rvs(sks.nb_r1[t][i],sks.nb_p1[t][i]) for i in range(nsg)] for t in range(nsample)]
  sks.nb_count=[[nbinom.rvs(sks.nb_r[t][i],sks.nb_p[t][i]) for i in range(nsg)] for t in range(nsample+1)]
  print('efficient: '+' '.join([str(x) for x in sks.isefficient]))
  return (sks)



def testcase1():
  '''
  Run test case 1
  '''
  sks=gentestcase1()
  iteratenb(sks)

def testcase2():
  '''
  Run test case 2
  '''
  sks=gentestcase2()
  iteratenb(sks)

def testcase3():
  '''
  Run test case 3
  '''
  sks=gentestcase3()
  iteratenbem(sks)

def wholegenometest1():
  '''
  Run the whole genome test
  '''
  # reading the file
  allgenedict=read_gene_from_file('/Users/wei/Dropbox/work/crispr/timdata/tim.norm.normalized.txt')
  desmat=np.matrix([[1,1,1,1],[0,0,1,0],[0,0,0,1]]).getT()
  for (gid, gene) in allgenedict.iteritems():
    gene.design_mat=desmat
  mycgene=allgenedict['MYC']
  iteratenbem(mycgene,estimateeff=True)
  return (allgenedict,mycgene)


def wholegenometest2(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  # reading the file
  pvargs='-k data/tim.norm.normalized.txt -d 1,0,0;1,0,0;1,1,0;1,0,1 -b "baseline,hl60,kbm7" -n results/tim/timtest'.split()
  
  rv=mageckmle_main(pvargs,returndict=True)
  return rv



def wholegenometest3(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  # reading the file
  allgenedict=read_gene_from_file('data/tim.norm.normalized.txt')
  desmat=np.matrix([[0,1,0],[0,0,1]]).getT()
  desmat=np.matrix([[1,1,1,1],[0,0,1,1],[0,0,0,1]]).getT()
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True)
    ngene+=1
    if ngene>maxgene:
      break
  write_gene_to_file(allgenedict,'results/tim/tim.kbm7beta.gene.txt')
  write_sgrna_to_file(allgenedict,'results/tim/tim.kbm7beta.sgrna.txt')
  return allgenedict


def wholegenometest4(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  # reading the file
  allgenedict=read_gene_from_file('data/tim.norm.normalized.txt')
  desmat=np.matrix([[0,1,0],[0,0,1]]).getT()
  desmat=np.matrix([[1,1,1,1],[0,0,1,1.34],[0,0,0,1.34]]).getT()
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True)
    ngene+=1
    if ngene>maxgene:
      break
  write_gene_to_file(allgenedict,'results/tim/tim.kbm7beta.1_34.gene.txt')
  write_sgrna_to_file(allgenedict,'results/tim/tim.kbm7beta.1_34.sgrna.txt')
  return allgenedict

def wholegenometest5(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  maxfittinggene=100
  # reading the file
  allgenedict=read_gene_from_file('data/tim.norm.normalized.txt')
  desmat=np.matrix([[1,1,1,1],[0,0,1,0],[0,0,0,1]]).getT()
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True)
    ngene+=1
    if ngene>maxfittinggene:
      break
  # model the mean and variance
  write_gene_to_file(allgenedict,'results/tim/tim.meanvar_initial.gene.txt')
  write_sgrna_to_file(allgenedict,'results/tim/tim.meanvar_initial.sgrna.txt')
  print('Modeling the mean and variance ...')
  mrm=MeanVarModel()
  mrm.get_mean_var_residule(allgenedict)
  mrm.model_mean_var_by_lm()
  # mrm.save_k_residule_to_file('results/tim/tim.meanvar.model.txt')
  
  # run the test again...
  print('Run the algorithm for the second time ...')
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True,meanvarmodel=mrm,restart=True)
    ngene+=1
    if ngene>maxgene:
      break
  # permutation
  iteratenbem_permutation(allgenedict,nround=100)
  # write to file
  write_gene_to_file(allgenedict,'results/tim/tim.meanvar.gene.txt')
  write_sgrna_to_file(allgenedict,'results/tim/tim.meanvar.sgrna.txt')
  return (allgenedict,mrm)


def wholegenometest_shalem(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  maxfittinggene=1000
  # reading the file
  allgenedict=read_gene_from_file('data/shalem.normalized.txt',includesamples=['plasmid','D7_R1','D7_R2','PLX7_R1','PLX7_R2'])
  desmat=np.matrix([[1,1,1,1,1],[0,1,1,0,0],[0,0,0,1,1]]).getT()
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True)
    ngene+=1
    if ngene>maxfittinggene:
      break
  # model the mean and variance
  write_gene_to_file(allgenedict,'results/shalem/shalem.meanvar_initial.gene.txt')
  write_sgrna_to_file(allgenedict,'results/shalem/shalem.meanvar_initial.sgrna.txt')
  print('Modeling the mean and variance ...')
  mrm=MeanVarModel()
  mrm.get_mean_var_residule(allgenedict)
  mrm.model_mean_var_by_lm()
  # mrm.save_k_residule_to_file('results/tim/tim.meanvar.model.txt')
  
  # run the test again...
  print('Run the algorithm for the second time ...')
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True,meanvarmodel=mrm,restart=True)
    ngene+=1
    if ngene>maxgene:
      break
  # write to file
  write_gene_to_file(allgenedict,'results/shalem/shalem.meanvar.gene.txt')
  write_sgrna_to_file(allgenedict,'results/shalem/shalem.meanvar.sgrna.txt')
  return (allgenedict,mrm)


def wholegenometest_shalem_plx7(maxgene=np.inf):
  '''
  Run the whole genome test
  '''
  maxfittinggene=1000
  # reading the file
  allgenedict=read_gene_from_file('data/shalem.normalized.txt',includesamples=['plasmid','D7_R1','D7_R2','PLX7_R1','PLX7_R2'])
  desmat=np.matrix([[1,1,1,1,1],[0,1,1,1,1],[0,0,0,1,1]]).getT()
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True)
    ngene+=1
    if ngene>maxfittinggene:
      break
  # model the mean and variance
  #write_gene_to_file(allgenedict,'results/shalem/shalem.meanvar_initial.gene.txt')
  #write_sgrna_to_file(allgenedict,'results/shalem/shalem.meanvar_initial.sgrna.txt')
  print('Modeling the mean and variance ...')
  mrm=MeanVarModel()
  mrm.get_mean_var_residule(allgenedict)
  mrm.model_mean_var_by_lm()
  # mrm.save_k_residule_to_file('results/tim/tim.meanvar.model.txt')
  
  # run the test again...
  print('Run the algorithm for the second time ...')
  ngene=0
  for (tgid,tginst) in allgenedict.iteritems():
    print('Calculating '+tgid+' ('+str(ngene)+') ... ')
    tginst.design_mat=desmat
    iteratenbem(tginst,debug=False,estimateeff=True,meanvarmodel=mrm,restart=True)
    ngene+=1
    if ngene>maxgene:
      break
  # write to file
  write_gene_to_file(allgenedict,'results/shalem/shalem.plx7beta.gene.txt')
  write_sgrna_to_file(allgenedict,'results/shalem/shalem.plx7beta.sgrna.txt')
  return (allgenedict,mrm)



def wholegenometest_sam(maxgene=np.inf):
  '''
  SAM test
  '''
  pvargs='-k data/sam_a375.normalized.txt -d 1,0,0;1,0,0;1,1,0;1,0,1;1,1,0;1,0,1  -n results/crispra/tmp.test -i zeo_plasmid_library,puro_plasmid_library,zeo_d3_rep_1,puro_d3_rep1,zeo_d3_rep2,puro_d3_rep2'.split()
  # Tracer()()
  rv=mageckmle_main(pvargs,returndict=True)
  return rv

def wholegenometest_t47d(maxgene=np.inf):
  '''
  SAM test
  '''
  pvargs='-k data/sam_a375.normalized.txt -d 1,0,0;1,0,0;1,1,0;1,0,1;1,1,0;1,0,1  -n results/crispra/tmp.test -i zeo_plasmid_library,puro_plasmid_library,zeo_d3_rep_1,puro_d3_rep1,zeo_d3_rep2,puro_d3_rep2'.split()
  pvargs='-k data/150131_seq.gecko.combined.txt -d results/t47d/designmatrix_t47d.txt -n results/t47d/t47d_all -i T47D_day0,T47D_w1_veh,T47D_w1_E2,T47D_w2_veh,T47D_w2_E2,T47D_w3_veh,T47D_w3_E2,T47D_w4_veh,T47D_w4_E2'.split()
  # Tracer()()
  rv=mageckmle_main(pvargs,returndict=True)
  return rv


def wholegenometest_tim(maxgene=np.inf):
  '''
  SAM test
  '''
  # pvargs='-k /home/wl948/datarun/wl948/project/crispr/mle/data_raw/tim.raw.data.new.txt -d /home/wl948/datarun/wl948/project/crispr/mle/designmat/designmat_tim.txt -n test --genes-var 0 --sgrna-efficiency /home/wl948/datarun/wl948/Dropbox/work/crispr/sgrnaeff/0.1/results/tim.library.ko.out --sgrna-eff-name-column 1 --sgrna-eff-score-column 3'.split()
  pvargs='-k /home/wl948/datarun/wl948/project/crispr/mle/data_raw/tim.raw.data.new.txt -d /home/wl948/datarun/wl948/project/crispr/mle/designmat/designmat_tim.txt -n test --genes-var 0'.split()
  # Tracer()()
  rv=mageckmle_main(pvargs,returndict=True)
  return rv





if __name__ == '__main__':
  try:
    #wholegenometest2();   
    #wholegenometest3();   
    # wholegenometest4();   
    #wholegenometest5(maxgene=100);   
    wholegenometest5()   
    #wholegenometest_shalem()
    #wholegenometest_shalem_plx7()
  except KeyboardInterrupt:
    sys.stderr.write("User interrupt me! ;-) Bye!\n")
    sys.exit(0)


