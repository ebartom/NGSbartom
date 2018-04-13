'''
MAGeCK MLE multiprocessing
'''

from __future__ import print_function
import re
import sys
import logging
import multiprocessing
import copy
import numpy as np

from mageck.mleem import iteratenbem

# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

def thread_p_func(dinst,args,iteratenbemargs,returndict):
  '''
  functions for multithreading
  Parameters:
    dist
        A dictionary of instances
    iteratenbemargs
        A dictionary of arguments in iteratenbem() function
        
  '''
  name = multiprocessing.current_process().name
  ngene=0
  logging.info(name+': total '+str(len(dinst))+ ' instances.')
  for (tgid,tginst) in dinst.iteritems():
    if ngene % 1000 ==1 or args.debug:
      logging.info(name+': Calculating '+tgid+' ('+str(ngene)+') ... ')
    iteratenbem(tginst,**iteratenbemargs)
    returndict[tgid]=tginst
    ngene+=1

  
def runem_multiproc(allgenedict,args,nproc=1, argsdict={}):
  '''
  Calling iternatembem using different number of threads
  Arguments:
    allgenedict:
        a dictionary of all gene instances
    args:
        arguments
    nproc
        The number of threads
    argsdict
        Positional arguments for iteratenbem
  '''
  # separate dicts
  instdictlist=[]
  mnger=multiprocessing.Manager()
  retdict=mnger.dict()
  if nproc==1:
    instdictlist.append(allgenedict)
  elif nproc<=0:
    logging.error('Error: incorrect number of threads.')
    sys.exit(-1)
  else:
    ngene=0
    instdictlist=[]
    for i in range(nproc):
      instdictlist.append({})
    for (tgid,tginst) in allgenedict.iteritems():
      targetlistid=ngene %nproc
      instdictlist[targetlistid][tgid]=tginst
      ngene+=1
  # start jobs
  jobs=[]
  
  for i in range(nproc):
    j=multiprocessing.Process(target=thread_p_func, name='Thread '+str(i),args=(instdictlist[i],args,argsdict,retdict))
    jobs.append(j)
    j.start()
    logging.info(j.name+' started.')
  
  for jj in jobs:
    jj.join()
    logging.info(jj.name+' completed.')
  logging.info('All threads completed.')
  # save the instance
  # Tracer()()
  for tgid in retdict.keys():
    tginst=retdict[tgid]
    allgenedict[tgid]=tginst


  
 

def iteratenbem_permutation(genedict,args,debug=True,nround=100,removeoutliers=False,size_factor=None):
  '''
  Perform permutation test
  '''
  logging.info('Start permuting '+str(nround)+' rounds ...')
  allsg=[]
  desmat=genedict[genedict.keys()[0]].design_mat
  nbeta1=desmat.shape[1]-1
  ngene=len(genedict)
  for (geneid, geneinst) in genedict.iteritems():
    nsg=geneinst.nb_count.shape[1]
    nsample=geneinst.nb_count.shape[0]
    countmat=geneinst.nb_count.getT()
    sgitem=[(geneinst.w_estimate[i],countmat[i]) for i in range(nsg)]
    allsg+=sgitem
  logging.info('Collecting '+str(len(allsg))+' sgRNAs from '+str(ngene)+' genes.')
  #
  genedictcopy=copy.deepcopy(genedict)
  betazeros=np.zeros((nround*ngene,nbeta1))
  #
  betaz_id=0
  for nrd in range(nround):
    np.random.shuffle(allsg)
    #
    logging.info('Permuting round '+str(nrd)+' ...')
    nid=0
    for (geneid, geneinst) in genedictcopy.iteritems():
      nsg=geneinst.nb_count.shape[1]
      nsample=geneinst.nb_count.shape[0]
      selitem=allsg[nid:nid+nsg]
      countmat=np.vstack([x[1] for x in selitem])
      w_es=np.array([x[0] for x in selitem])
      geneinst.nb_count=countmat.getT()
      geneinst.w_estimate=w_es
      nid+=nsg
    # end gene loop
    #iteratenbem(geneinst,debug=False,estimateeff=True,updateeff=False,removeoutliers=removeoutliers,size_factor=size_factor,logem=False)
    argsdict={'debug':False,'estimateeff':True,'updateeff':False,'removeoutliers':removeoutliers,'size_factor':size_factor,'logem':False}
    runem_multiproc(genedictcopy,args,nproc=args.threads,argsdict=argsdict)
    for (geneid, geneinst) in genedictcopy.iteritems():
      nsg=geneinst.nb_count.shape[1]
      beta_es=geneinst.beta_estimate[nsg:]
      betazeros[betaz_id,:]=beta_es
      betaz_id+=1
    # end gene loop
  # end permutation
  logging.info('Assigning p values...')
  ncompare=betazeros.shape[0]*1.0
  for (geneid, geneinst) in genedict.iteritems():
    nsg=geneinst.nb_count.shape[1]
    beta_es=geneinst.beta_estimate[nsg:]
    cp_u0=np.sum(betazeros>beta_es,axis=0)
    cp_u1=np.sum(betazeros<beta_es,axis=0)
    cp_ustack=np.vstack((cp_u0/ncompare,cp_u1/ncompare))
    cp_minval=np.min(cp_ustack,axis=0)
    #cp_minvec=np.array(cp_minval)[0]
    cp_minvec=cp_minval*2
    geneinst.beta_permute_pval=cp_minvec
    geneinst.beta_permute_pval_neg=cp_ustack[1]
    geneinst.beta_permute_pval_pos=cp_ustack[0]
    # Tracer()()

  return betazeros


