'''
Defining the core EM MLE approach
'''

from __future__ import print_function

import re
import sys
import scipy
from scipy.stats import nbinom
from scipy.stats import norm
import random
import math
import numpy as np
import numpy.linalg as linalg
import copy

from mageck.mleclassdef import *
from mageck.mledesignmat import *
from mageck.mlemeanvar import *
# from mageck.mlemultiprocessing import runem_multiproc

import logging


# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

def getloglikelihood(kmat,mu_estimate,alpha):
  '''
  Get the log likelihood estimation of NB, using the current estimation of beta
  '''
  # logmu_est=sk.extended_design_mat * np.matrix(beta_est).getT()
  # these are all N*1 matrix
  mu_vec=[t[0] for t in mu_estimate.tolist()]
  k_vec=[t[0] for t in kmat.tolist()]
  if len(mu_vec) != len(k_vec):
    raise ValueError('Count table dimension is not the same as mu vector dimension.')
  var_vec=[t+alpha*t*t for t in mu_vec]
  nb_p=[mu_vec[i]/var_vec[i] for i in range(len(mu_vec))] 
  nb_r=[mu_vec[i]*mu_vec[i]/(var_vec[i]-mu_vec[i]) for i in range(len(mu_vec))]
  logp=[nbinom.logpmf(k_vec[i],nb_r[i],nb_p[i]) for i in range(len(mu_vec))]
  return sum(logp)

def getloglikelihood2(kmat,mu_estimate,alpha,sumup=False,log=True):
  '''
  Get the log likelihood estimation of NB, using the current estimation of beta
  '''
  #logmu_est=sk.extended_design_mat * np.matrix(beta_est).getT()
  # Tracer()()
  #mu_estimate= np.exp(logmu_est)
  # these are all N*1 matrix
  #mu_vec=np.array([t[0] for t in mu_estimate.tolist()])
  #k_vec=np.array([round(t[0]) for t in kmat.tolist()])
  #if len(mu_vec) != len(k_vec):
  #  raise ValueError('Count table dimension is not the same as mu vector dimension.')
  # var_vec=mu_vec+alpha*mu_vec*mu_vec
  # nb_p=[mu_vec[i]/var_vec[i] for i in range(len(mu_vec))] 
  # nb_r=[mu_vec[i]*mu_vec[i]/(var_vec[i]-mu_vec[i]) for i in range(len(mu_vec))]
  # if log:
  #  logp=np.array([nbinom.logpmf(k_vec[i],nb_r[i],nb_p[i]) for i in range(len(mu_vec))])
  #else:
  #  logp=np.array([nbinom.pmf(k_vec[i],nb_r[i],nb_p[i]) for i in range(len(mu_vec))])
  if kmat.shape[0] != mu_estimate.shape[0]:
    raise ValueError('Count table dimension is not the same as mu vector dimension.')
  kmat_r=np.round(kmat)
  mu_sq=np.multiply(mu_estimate,mu_estimate)
  var_vec=mu_estimate+np.multiply(alpha, mu_sq)
  nb_p=np.divide(mu_estimate,var_vec)
  nb_r=np.divide(mu_sq,var_vec-mu_estimate)
  if log:
    logp=nbinom.logpmf(kmat_r,nb_r,nb_p)
  else:
    logp=nbinom.pmf(kmat,nb_r,nb_p)
  
  if np.isnan(np.sum(logp)):
    #raise ValueError('nan values for log likelihood!')
    logp=np.where(np.isnan(logp),0,logp)
  if sumup:
    return np.sum(logp)
  else:
    return logp

def remove_outlier(y_residule,w_matrix,outprob=0.05):
  '''
  During ridge regression, remove outliers.
  Diagnoal values of w_matrix whose corresponding y_residules values above the outprob (%) and below the outprob (%) are set to 0.
  Currently, we only remove the maximum and minimum values; but in the future, all values with max or min outprob (%) should be set to 0.
  '''
  ny=y_residule.shape[0]
  nsel=int(ny*outprob)
  if nsel<1:
    return 0
  threshold=2.0
  # choose the top nsel elements
  npmaxi=np.argmax(y_residule)
  npmini=np.argmin(y_residule)
  nmax_val=np.sum(y_residule>y_residule[npmaxi]/threshold)
  if nmax_val==1:
    w_matrix[npmaxi,npmaxi]=0
  nmin_val=np.sum(y_residule<y_residule[npmini]/threshold)
  if nmin_val==1:
    w_matrix[npmini,npmini]=0
  return 0
  
def em_whileloop(sk,beta_init_mat_0,size_vec,wfrac_list_0,alpha_dispersion,alpha_val,estimateeff,updateeff,removeoutliers,debug):
  '''
  The while loop in EM
  Required variables:
    beta_init_mat_0: beta_init_mat before calling em_whileloop
    size_vec
    wfrac_list_0: wfrac_list before calling em_wileloop
    alpha_dispersion
    alpha_val (temporary)
    estimateeff
    updateeff
    removeoutliers
    debug

    sk.extdesign_mat
    sk.datak_mat
    sk.nb_count
  Return value:
    ret: an integer to indicate the status: 
      0  success (converges)
      1  success (but it reaches the max number of iterations)
      2  success (has a nan value on beta, so returns to a last iterated beta that does not equal to nan; an indication that there are possible numeric problems)
    Other variables:
      wfrac_list
      beta_new
      beta1_se_mat
      mu_estimate
      sgrna_residule
  '''
  # constants
  n_max_init=1000
  diff_cutoff=1e-9
  max_abs_beta=100 # the maximum abs value of beta
  # preparing variables for loop
  beta_init_mat=beta_init_mat_0.copy()
  wfrac_list=wfrac_list_0.copy()
  #extdesign_mat=sk.extended_design_mat
  #extdesignmat_residule=sk.extended_design_mat_residule
  datak_mat=sk.sgrna_kvalue
  n=(sk.nb_count.shape[1]) # the number of sgRNAs
  nallsample=sk.design_mat.shape[0]
  nsample=sk.design_mat.shape[0]-1 # the number of samples excluding the 1st base sample
  nbeta1=sk.design_mat.shape[1]-1 # the number of betas excluding the 1st beta (baseline beta)
  (bsid,design_mat,extdesign_mat,extdesignmat_residule)=DesignMatCache.get_record(n)
  if np.isscalar(alpha_dispersion):
    alpha_dispersion_mat=alpha_dispersion
  else:
    alpha_dispersion_mat=np.matrix(alpha_dispersion).T
  n_iter=1
  #
  retval=0
  #
  # intermediate variables for calculating wald test statistics
  beta1_se_mat=None
  beta1_new=None
  # while loop
  while(True):
    logmu_estimate=extdesign_mat*beta_init_mat # (nsample*nsgrna)*1 matrix
    mu_estimate= np.multiply(size_vec,np.exp(logmu_estimate))
    sgrna_residule=(datak_mat-mu_estimate)
    if estimateeff and True:
      w_coef=np.matrix([[1.0]]*n+[[x] for x in wfrac_list]*nsample+[[1.0-x] for x in wfrac_list]*nsample)
      z_estimate=np.multiply(sgrna_residule,w_coef)/mu_estimate + logmu_estimate
    else:
      z_estimate=sgrna_residule/mu_estimate + logmu_estimate
    # calculate log likelihood
    # Tracer()()
    logl=getloglikelihood2(datak_mat,mu_estimate,alpha_dispersion_mat)
    if np.isnan(np.sum(logl)):
      raise ValueError('nan values for log likelihood!')
    #nologl=getloglikelihood2(datak_mat,sk,beta_init_mat.getT(),alpha_val,log=False)
    logll=np.sum(logl) #logll_list+=[logll]
    # update posterior prob. that an sgRNA is efficient
    if estimateeff and updateeff:
      # iteratively update efficiency
      beta_0_mat_ext=beta_init_mat.copy()
      beta_0_mat_ext[n:n+nbeta1,:]=0
      logmu_estimate_0=extdesign_mat*beta_0_mat_ext # (nsample*nsgrna)*1 matrix
      mu_estimate_0= np.multiply(size_vec,np.exp(logmu_estimate_0))
      logl_0=getloglikelihood2(datak_mat,mu_estimate_0,alpha_dispersion_mat)
      #nologl_0=getloglikelihood2(datak_mat,sk,beta_0_mat_ext.getT(),alpha_val,log=False)
      #logl_0=getloglikelihood2(datak_mat,sk,beta_0_mat_ext.getT(),alpha_dispersion)
      # use traditional no-log p value method
      #nolog_sum=nologl+nologl_0
      #for i in range(nolog_sum.shape[0]):
      #  if nolog_sum[i]==0:
      #    nolog_sum[i]=1
      #sgfrac=nologl/nolog_sum
      #sgfrac_mat=np.matrix(sgfrac)
      # use log p value
      if False:
        sgfrac_mat=np.matrix(1/(1+np.exp(logl_0-logl)))
        # sgfrac_mat_res=sgfrac_mat.reshape(2*nsample+1,n).getT()
        # sgfrac_mat_res=sgfrac_mat[0,0:(n*(nsample+1))].reshape(1*nsample+1,n).getT()
        sgfrac_mat_res=sgfrac_mat[n:(n*(nsample+1)),:].reshape(1*nsample,n).getT()
        wfrac=np.mean(sgfrac_mat_res,axis=1)
      else:
        sgfrac_mat=np.matrix(logl_0-logl)
        sgfrac_mat_res=sgfrac_mat[n:(n*(nsample+1)),:].reshape(1*nsample,n).getT()
        sgfrac_mat_res=np.where(sgfrac_mat_res>100,100,sgfrac_mat_res)
        ## here, we should use a minimum value to estimate the efficiency, instead of the mean value. In another way, if one sgRNA is shown to be efficient in only 1 conditions, it should be efficient
        # wfrac_diff=np.mean(sgfrac_mat_res,axis=1)
        # Tracer()()
        wfrac_diff=np.min(sgfrac_mat_res,axis=1)
        #
        # wfrac_diff=np.where(wfrac_diff>200,200,wfrac_diff)
        wfrac=np.matrix(1/(1+np.exp(wfrac_diff)))
      wfrac_list=wfrac.getA1()
      if debug:
        print('frac:'+' '.join([decformat(x) for x in wfrac.getT().tolist()[0]]))
    else:
      # don't estimate sgrna efficiency
      # wfrac=np.matrix([[1]])
      pass
    # end if  
    
    # w
    #w_list=np.array([t[0] for t in mu_estimate.tolist()]); # dimension: (nsg*nsample)*1 
    w_list=mu_estimate.getA1() # dimension: (nsg*nsample)*1 
    # Tracer()()
    if estimateeff:
      # modify the w
      wfrac_ext=np.append(np.tile(wfrac_list,(nsample+1)),np.tile(1-wfrac_list,nsample))
      #wfrac_ext_bin=wfrac_ext>=1e-4
      #wfrac_ext=wfrac_ext*wfrac_ext_bin+(1e-4)*(1-wfrac_ext_bin)
      wfrac_ext=np.where(wfrac_ext>=1e-2,wfrac_ext,1e-2)
      #for i in range(len(wfrac_ext)):
      #  if wfrac_ext[i]<1e-4:
      #    wfrac_ext[i]=1e-4
      #w_list_ele=[1.0/(1.0/w_list[ti]+alpha_val) for ti in range(len(w_list))]
      #for i in range(n,len(w_list_ele)):
      #  w_list_ele[i]=1.0*wfrac_ext[i]/(1.0/w_list[i]+alpha_val)
      w_list_ele=1.0/(1.0/w_list+alpha_dispersion)
      w_list_ele[n:]=w_list_ele[n:]*wfrac_ext[n:]
      # do not modify the w
      # w_list_ele=[1/(1/t+alpha_val) for t in w_list]
    else:
      # w_list_ele=[1/(1/t+alpha_val) for t in w_list]
      w_list_ele=1.0/(1.0/w_list+alpha_dispersion)
    # the following heuristic step makes sure the w value doesn't go too small (and beta values don't go too large)
    # Tracer()()
    wlmax=np.max(w_list_ele)/100
    w_list_ele=np.where(w_list_ele<wlmax,wlmax,w_list_ele)
    if debug:
        print('w:'+' '.join([str(x) for x in w_list_ele]))
    w_matrix=np.diag(w_list_ele)
    #
    # inv 
    if True:
      # solution 0: directly solve z_estimate=extdesign_mat*beta
      xwx_mat=extdesign_mat.getT()*w_matrix*extdesign_mat+alpha_val*np.matrix(np.identity(n+nbeta1))
      xwx_inv=linalg.inv(xwx_mat)
      beta_new=xwx_inv*extdesign_mat.getT()*w_matrix*z_estimate
      # matrix for Standard Error
      beta1_new=beta_new[n:,:]
      beta1_new_pref=xwx_inv*extdesign_mat.getT()*w_matrix
      beta1_se_mat=beta1_new_pref*extdesign_mat*xwx_inv
      beta1_se_mat=beta1_se_mat[n:,n:]
    if False:
      # solution 1: uses the iteratively updated least squares 
      # update both beta0 and beta1
      z_residule=z_estimate-extdesign_mat[:,0:n]*beta_init_mat[0:n,:]
      #z_residule=z_residule[0:(nsample+1)*n,:] # only consider samples from efficient sgRNAs
      #w_matrix_residule=w_matrix[0:(nsample+1)*n,0:(nsample+1)*n]
      #
      # remove outliers: set up the w_matrix
      if removeoutliers:
        remove_outlier(z_residule,w_matrix)
      # perform the matrix operation
      xwx_mat=extdesign_mat.getT()*w_matrix*extdesign_mat+alpha_val*np.matrix(np.identity(extdesign_mat.shape[1]))
      xwx_inv=linalg.inv(xwx_mat)
      beta1_new_pref=xwx_inv*extdesign_mat.getT()*w_matrix
      beta1_new=beta1_new_pref*z_residule
      beta_new=beta_init_mat.copy()
      # Tracer()()
      
      beta_new=beta1_new
      #beta_new[n:,:]=beta1_new
      # matrix for Standard Error
      beta1_se_mat=beta1_new_pref*extdesign_mat*xwx_inv
    if False:
      # solution 2: uses the iteratively updated least squares 
      # keep the beta0, only update beta1
      z_residule=z_estimate-extdesign_mat[:,0:n]*beta_init_mat[0:n,:]
      z_residule=z_residule[0:(nsample+1)*n,:] # only consider samples from efficient sgRNAs
      w_matrix_residule=w_matrix[0:(nsample+1)*n,0:(nsample+1)*n]
      #
      # remove outliers: set up the w_matrix
      if removeoutliers:
        remove_outlier(z_residule,w_matrix)
      # perform the matrix operation
      xwx_mat=extdesignmat_residule.getT()*w_matrix_residule*extdesignmat_residule+alpha_val*np.matrix(np.identity(design_mat.shape[1]))
      xwx_inv=linalg.inv(xwx_mat)
      beta1_new_pref=xwx_inv*extdesignmat_residule.getT()*w_matrix_residule
      beta1_new=beta1_new_pref*z_residule
      beta_new=beta_init_mat.copy()
      # Tracer()()
      
      beta_new[n:,:]=beta1_new
      # matrix for Standard Error
      beta1_se_mat=beta1_new_pref*extdesignmat_residule*xwx_inv
      ## Wald test p values
      ## moved to the end of the calculation
      # beta_se_val=(np.diag(beta1_se_mat))
      # beta_new_zscore=beta1_new.tolist()/beta_se_val
      # beta_new_zscore_list=[x[0] for x in beta_new_zscore.tolist()]
      # beta_new_pval=norm.sf(beta_new_zscore_list)*2
    #
    beta_diff=beta_new-beta_init_mat
    if np.isnan(np.sum(beta_new)) or np.max(np.abs(beta_new))>max_abs_beta:
      retval=2 # encounter some numeric problems
      break
    else:
      beta_init_mat=beta_new
    # logll=getloglikelihood2(datak_mat,sk,beta_init_mat.getT(),alpha_val)
    n_iter+=1
    ## if need break
    #diffval=sum([t[0]*t[0] for t in beta_diff.tolist()[n:]])
    #absval=sum([t[0]*t[0] for t in beta_new.tolist()[n:]])
    diffval=np.sum(beta_diff.getA1()[n:]**2)
    absval=np.sum(beta_new.getA1()[n:]**2)
    #Tracer()()
    if abs(absval)<1e-9:
      absval=1.0
    difffrac=diffval/absval
    # print information
    if debug:
      print('Iteration '+str(n_iter)+', updated beta:',end='')
      print(' '.join([decformat(x[0]) for x in beta_new.tolist()]))
      print('log likelihood:'+str(logll)+', frac:'+decformat(difffrac)+',abs:'+decformat(absval))
    if difffrac<diff_cutoff: 
      break
    if n_iter>n_max_init:
      retval=1
      break
  # end while
  # construct return value
  beta_return=beta_init_mat # equivalent to beta_new, but equals to last iteration beta value when beta_new=nan
  return (retval,wfrac_list,beta_return,beta1_se_mat,mu_estimate,sgrna_residule)


def iteratenbem(sk,debug=True,estimateeff=False,updateeff=True,plot=False,alpha_val=0.01,meanvarmodel=None,restart=False,removeoutliers=False,size_factor=None,logem=True):
  '''
  Iteratively solve the value of betas using negative binomial joint likelihood function.
  should work if estimateeff is false; for EM, use a new function.
  
  Parameters
  ----------
  sk
    The SimCaseSimple class of genes
  debug
    Whether to print debug information
  estimateeff
    Whether to use sgRNA efficiency. If set to false, don't use efficiency (i.e., all sgRNAs are efficient). If set to true: (1) if the value of sk.w_estimate is empty, assume all sgRNAs are efficient; (2) otherwise, use the values from sk.w_estimate as the initial value.
  updateeff
    Whether to iteratively update sgRNA efficiency. If False, always use the initial estimation of sk.w_estimate
  plot
    Whether to plot the diagnosis information for each EM iteration
  alpha_val
    Prior for beta
  meanvarmodel
    Model for Mean and Variance; default is None
  restart
    Whether the beta value should come from the previous iterations of EM
  removeoutliers
    Whether to remove outliers
  size_factor
    Size factor correcting for sequencing depth; default is None (assuming all sizes are 1)
  logem
    Whether to log detailed parameter adjustments for EM
  '''
  # parameters
  max_dim=10000 # the maximum number of variables defined to avoid memory error
  #n_max_init=50
  # constants
  n=(sk.nb_count.shape[1]) # the number of sgRNAs
  nallsample=sk.design_mat.shape[0]
  nsample=sk.design_mat.shape[0]-1 # the number of samples excluding the 1st base sample
  nbeta1=sk.design_mat.shape[1]-1 # the number of betas excluding the 1st beta (baseline beta)
  logll_list=[] # log likelihood
  
  if n*nsample>max_dim:
    # to avoid memory error, skip the instance
    sk.beta_estimate=np.array([0.0]*(nbeta1+n))
    sk.beta_pval=np.array([1.0]*(nbeta1))
    sk.beta_zscore=np.array([0.0]*(nbeta1))
    sk.beta_pval_pos=np.array([1.0]*(nbeta1))
    sk.beta_pval_neg=np.array([1.0]*(nbeta1))
    return

  # design mat
  # calculate the number of base line samples
  #(basesampleid,_new_designmat)=analyze_designmat(sk.design_mat)
  #design_mat=sk.design_mat[1:,1:]; # assuming the 1st column must be the baseline condition, and the 1st row is the base line condition; they are removed.
  #extdesign_mat=getextenddesignmat(n,nsample,design_mat,includebase=True)
  #extdesignmat_residule=extdesign_mat[0:(nsample+1)*n,n:]; # this is the design matrix containing only non-baseline betas and first (nsample+1) samples
  if DesignMatCache.has_record(n) == False:
    DesignMatCache.save_record(sk.design_mat,n)
  if DesignMatCache.has_record(n) == True:
    (basesampleid,design_mat,extdesign_mat,extdesignmat_residule)=DesignMatCache.get_record(n)
  else:
    raise ValueError('There is no corresponding record in DesignMatCache.')
  
  sk.extended_design_mat=extdesign_mat
  sk.extended_design_mat_residule=extdesignmat_residule
  
  if size_factor==None:
    size_vec=np.ones(n*(2*nallsample-1))
  else:
    if len(size_factor) != nallsample:
      raise ValueError('The provided size factor length does not equal to the number of samples.')
    size_vec=np.repeat(np.append(size_factor,size_factor[1:]),n)
  size_mat=np.matrix(size_vec.reshape([(2*nallsample-1),n]))
  size_mat=size_mat[:nallsample,:]
  size_vec=np.matrix(size_vec).getT()

  # prepare for the data matrix
  # old method: slow 
  #for i in range(nsample+1):
  #  dk_list+=sk.nb_count[i]
  #for i in range(1,nsample+1): # repeat the treatment sample sequences
  #  dk_list+=sk.nb_count[i]
  #datak=np.array(dk_list)
  #datak_mat=np.matrix(datak).getT()
  # convert to matrix operation
  datak_mat_0=np.vstack((sk.nb_count,sk.nb_count[1:,:]))
  datak_mat=datak_mat_0.reshape(datak_mat_0.shape[0]*datak_mat_0.shape[1],1)
  datak=datak_mat
  sk.sgrna_kvalue=datak_mat



  # prepare for the dispersion estimation
  if meanvarmodel == None:
    alpha_dispersion=alpha_val
    alpha_dispersion_mat=alpha_val
  else:
    # lm method
    alpha_dispersion_mat=np.matrix(meanvarmodel.get_lm_var(datak,returnalpha=True))
    # glm method
    #normalized_k=[]
    #inverse_size_f=[1/i for i in size_factor]
    #sg_k=[x[0] for x in sk.sgrna_kvalue.tolist()]
    ##logging.info(sg_k)
    #for i in range(len(sk.sgrnaid)):
    #    #logging.info(i)
    #    #logging.info([sg_k[i+k*n] for k in range(nallsample)])
    #    normalized_k_mean=np.mean(np.multiply(inverse_size_f,[sg_k[i+k*n] for k in range(nallsample)]))
    #    #normalized_k_mean=np.mean(np.multiply(inverse_size_f,sg_k[i*nallsample:(i+1)*nallsample]))
    #    normalized_k.append([normalized_k_mean])
    #normalized_k=normalized_k*nallsample
    #Tracer()
    #normalized_datak=np.matrix(normalized_k).reshape(datak_mat_0.shape[0]*datak_mat_0.shape[1],1)
    #alpha_dispersion_mat=np.matrix(meanvarmodel.get_glm_dispersion(normalized_datak,returnalpha=True))
    #alpha_dispersion_mat=np.matrix(meanvarmodel.get_glm_dispersion(datak,returnalpha=True))
    alpha_dispersion=alpha_dispersion_mat.getA1()
  # initial guess of beta
  baseline_sample_matrix=np.log(sk.nb_count[basesampleid,:])-np.log(size_mat[basesampleid,:]) 
  if restart==False or len(sk.beta_estimate)==0:
    beta_vec1=np.mean(baseline_sample_matrix,axis=0)
    # beta_vec1=beta_vec1.tolist()[0]
    # beta_vec1=[math.log(sk.nb_count[0][i]) for i in range(n)]
    beta_vec2=np.log(sk.nb_count[1:,:])-np.log(size_mat[1:,:])-beta_vec1
    #beta_vec2=np.array([[math.log(sk.nb_count[t][i]) - beta_vec1[i]  for i in range(n)] for t in range(1,nsample+1)])
    beta1_meanval=np.mean(beta_vec2,axis=1)
    # beta1_meanval=np.array([sum(beta_vec2[t])/n for t in range(nsample)])
    ## solving bmean = designmat * beta_es 
    beta1_es_mat=linalg.inv(design_mat.getT()*design_mat+alpha_val*np.matrix(np.identity(design_mat.shape[1])))*design_mat.getT()*(beta1_meanval)
    # beta1_es=np.array([t[0] for t in beta1_es_mat.tolist()])
    # only estimate beta1_es , not the whole beta
    # beta_init_mat=np.matrix(beta1_es).getT()
    # beta_init_mat=np.matrix(np.append(beta_vec1,beta1_es)).getT()
    beta_init_mat=np.vstack((beta_vec1.getT(),beta1_es_mat))
    if len(sk.w_estimate)==0 or estimateeff == False:
      # no prior information of w_estimate
      # wfrac_list=np.array([1.0]*n); # a list of size nsgRNA 
      wfrac_list=np.ones(n) # a list of size nsgRNA 
    else:
      # w_estimate was set up previously. Continue using them
      wfrac_list=np.array(sk.w_estimate)
  else:
    # only work if restart=True
    beta_init_mat=np.matrix([[x] for x in sk.beta_estimate])
    wfrac_list=np.array(sk.w_estimate)
  wfrac=0.0
  
  # initial log likelihood
  logmu_estimate=extdesign_mat*beta_init_mat # (nsample*nsgrna)*1 matrix
  mu_estimate= np.multiply(size_vec,np.exp(logmu_estimate))
  # Tracer()()
  logl=getloglikelihood2(datak_mat,mu_estimate,alpha_dispersion_mat)
  logll=sum(logl) 
  logll_list+=[logll]
  if debug:
    print('Initial log likelihood:'+str(logll))
    print('Initial beta:'+' '.join([decformat(x[0]) for x in beta_init_mat.tolist()]))
  
  #
  # start while loop
  ntrial=0
  while True:
    retlist=em_whileloop(sk,beta_init_mat,size_vec,wfrac_list,alpha_dispersion,alpha_val,estimateeff,updateeff,removeoutliers,debug)
    retval=retlist[0]
    wfrac_list=retlist[1]
    beta_new=retlist[2]
    beta1_se_mat=retlist[3]
    mu_estimate=retlist[4]
    sgrna_residule=retlist[5]
    if retval == 0:
      break
    elif retval == 1:
      if logem:
        logging.warning(sk.prefix+': reaches the maximum number of iterations.')
      #break
    elif retval == 2:
      if logem:
        logging.warning(sk.prefix+': beta value does not converge. Try to increase the value of alpha ..')
    if logem:
      logging.warning(sk.prefix+': alpha: '+str(alpha_val))
      #logging.warning(sk.prefix+': alpha_dispersion: '+str(alpha_dispersion))
    alpha_dispersion*=3
    alpha_val*=3
    ntrial+=1
    if ntrial > 10:
      if logem:
        logging.warning(sk.prefix+': after 10 tries, the beta value still does not converge. Exit ...')
      break
  # end while
  # save the result
  # sk.w_estimate=np.array([x for x in wfrac.getT().tolist()[0]])
  # sk.beta_estimate=np.array([x[0] for x in beta_init_mat.tolist()])
  if updateeff:
    sk.w_estimate=wfrac_list
  sk.beta_estimate=beta_new.getA1()
  # calculate p value
  beta_se_val=np.sqrt(np.diag(beta1_se_mat))
  #beta_new_zscore=beta1_new.tolist()/beta_se_val
  #beta_new_zscore_list=np.array([x[0] for x in beta_new_zscore.tolist()])
  beta1_new=beta_new[n:,:]
  beta_new_zscore_list=beta1_new.getA1()/beta_se_val
  # Tracer()()
  beta_new_pval=norm.sf(np.abs(beta_new_zscore_list))*2
  sk.beta_pval=beta_new_pval
  sk.beta_zscore=beta_new_zscore_list
  sk.beta_pval_pos=norm.sf(beta_new_zscore_list)
  sk.beta_pval_neg=norm.cdf(beta_new_zscore_list)
  sk.mu_estimate=mu_estimate
  sk.sgrna_residule=sgrna_residule
  if debug:
    print('Real beta:'+' '.join([decformat(x) for x in sk.beta0+sk.beta1]))
    print('z score:'+' '.join([decformat(x) for x in beta_new_zscore_list]))
    print('p value:'+' '.join([decformat(x) for x in beta_new_pval]))
  if plot:
    betaval=[sk.beta0+sk.beta1,sk.beta_estimate]
    if estimateeff:
      effvalue=[sk.isefficient]
      effvalue+=[sk.w_estimate]
      plotem(logll_list,betaval,effval=effvalue,filename=sk.prefix)
    else:
      plotem(logll_list,betaval,filename=sk.prefix)
   
def plotem(logll,betaval,effval=None,filename='sample1'):
  '''
  Plot figures for EM
  '''
  import matplotlib.pyplot as plt
  from matplotlib.backends.backend_pdf import PdfPages
  pp = PdfPages(filename+'.pdf')
  
  fig1=plt.figure(1)
  plt.subplot(211)
  plt.plot(logll,'bo-')
  plt.xlabel('Iterations')
  plt.ylabel('log likelihood')

  plt.subplot(212)
  index=np.arange(len(betaval[0]))
  bwidth=0.35
  realval=plt.bar(index,betaval[0],bwidth,color='b',label='True')
  simvar=plt.bar(index+bwidth,betaval[1],bwidth,color='r',label='Estimated')
  plt.legend()

  plt.xlabel(r'$\beta$')
  plt.xticks(index + bwidth, [str(x+1) for x in index])
  # plt.show()
  pp.savefig(fig1)
  
  if effval != None:
    fig2=plt.figure(2)
    index=np.arange(len(effval[0]))
    bwidth=0.35
    realval=plt.bar(index,effval[0],bwidth,color='b',label='True')
    simvar=plt.bar(index+bwidth,effval[1],bwidth,color='r',label='Estimated')
    plt.legend()
    plt.xlabel('Efficiency estimation')
    plt.xticks(index + bwidth, [str(x+1) for x in index])
    plt.axis([0,len(effval[0]),0,1.2])
    pp.savefig(fig2)
  # plt.show()
  pp.close()
 


