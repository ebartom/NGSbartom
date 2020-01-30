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

from mleclassdef import *
from mledesignmat import *
import logging

def getloglikelihood2(kmat,mu_estimate,alpha,sumup=False,log=True):
    '''
    Get the log likelihood estimation of NB, using the current estimation of beta
    '''
    if kmat.shape[0] != mu_estimate.shape[0]:
        raise ValueError('Count table dimension is not the same as mu vector dimension.')
    alpha=np.matrix(alpha).reshape(mu_estimate.shape[0],mu_estimate.shape[1])
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

def getloglikelihood3(kmat,mu_estimate,alpha,sumup=False,log=True):
    '''
    Get the log likelihood estimation of NB, using the current estimation of beta
    '''
    if kmat.shape[0] != mu_estimate.shape[0]:
        raise ValueError('Count table dimension is not the same as mu vector dimension.')
    alpha=np.matrix(alpha).reshape(mu_estimate.shape[0],mu_estimate.shape[1])
    kmat_r=np.round(kmat)
    mu_sq=np.multiply(mu_estimate,mu_estimate)
    var_vec=mu_estimate+np.multiply(alpha, mu_sq)
    nb_p=np.divide(mu_estimate,var_vec)
    nb_r=np.divide(mu_sq,var_vec-mu_estimate)
    p=nbinom.cdf(kmat,nb_r,nb_p)
    p=np.where(p<0.5,p,1-p)
    if log:
        #logp=nbinom.logcdf(kmat_r,nb_r,nb_p)
        logp=np.log(p)
    else:
        logp=p
        #logp=nbinom.cdf(kmat,nb_r,nb_p)
#

    if np.isnan(np.sum(logp)):
        #raise ValueError('nan values for log likelihood!')
        logp=np.where(np.isnan(logp),0,logp)
    if sumup:
        return np.sum(logp)
    else:
        return logp

def iteratenbem(sk,debug=True,estimateeff=False,updateeff=True,plot=False,PPI_prior=False,alpha_val=0.01,meanvarmodel=None,restart=True,removeoutliers=False,size_factor=None,beta1_prior_var=None):
    # parameters
    n_max_init=1000
    diff_cutoff=1e-9

    n=(sk.nb_count.shape[1])
    nallsample=sk.design_mat.shape[0]
    nsample=sk.design_mat.shape[0]-1 # the number of samples excluding the 1st base sample
    nbeta1=sk.design_mat.shape[1]-1 # the number of betas excluding the 1st beta (baseline beta)
    logll_list=[] # log likelihood
    if DesignMatCache.has_record(n) == False:
        DesignMatCache.save_record(sk.design_mat,n)
    if DesignMatCache.has_record(n) == True:
        (basesampleid,design_mat,extdesign_mat,extdesignmat_residule)=DesignMatCache.get_record(n)
    else:
        raise ValueError('There is no corresponding record in DesignMatCache.')

    extdesign_mat=extdesign_mat[:(sk.nb_count.shape[0])*(sk.nb_count.shape[1]),]
    sk.extended_design_mat=extdesign_mat
    #logging.info(sk.nb_count.shape)
    #logging.info(extdesign_mat)
    if PPI_prior==False:
        if beta1_prior_var!=None:
            beta1_prior_inverse_variance=[1/i for i in beta1_prior_var]
        else:
            beta1_prior_inverse_variance=[0]*nbeta1
    else:
        beta1_prior_inverse_variance=[1/i for i in sk.prior_variance]
        beta_prior_mean=np.matrix([0]*n+sk.prior_mean)
        beta_prior_mean=beta_prior_mean.reshape(n+len(sk.prior_mean),1)
    beta_prior_inverse_variance=[0]*n+beta1_prior_inverse_variance

    if size_factor==None:
        size_vec=np.ones(n*(nallsample))
    else:
        if len(size_factor) != nallsample:
            raise ValueError('The provided size factor length does not equal to the number of samples.')
        size_vec=np.repeat(size_factor,n)

    size_mat=np.matrix(size_vec.reshape([nallsample,n]))
    size_vec=np.matrix(size_vec).getT()

    datak_mat_0=sk.nb_count
    #logging.info(datak_mat_0)
    datak=datak_mat_0.reshape(datak_mat_0.shape[0]*datak_mat_0.shape[1],1)
    sk.sgrna_kvalue=datak
    #logging.info(datak)
    normalized_k=[]
    inverse_size_f=[1/i for i in size_factor]
    sg_k=[x[0] for x in sk.sgrna_kvalue.tolist()]
    #logging.info(sg_k)
    for i in range(len(sk.sgrnaid)):
        #logging.info(i)
        #logging.info([sg_k[i+k*n] for k in range(nallsample)])
        normalized_k_mean=np.mean(np.multiply(inverse_size_f,[sg_k[i+k*n] for k in range(nallsample)]))
        #normalized_k_mean=np.mean(np.multiply(inverse_size_f,sg_k[i*nallsample:(i+1)*nallsample]))
        normalized_k.append([normalized_k_mean])
    normalized_k=normalized_k*nallsample
    normalized_datak=np.matrix(normalized_k).reshape(datak_mat_0.shape[0]*datak_mat_0.shape[1],1)

    # prepare for the dispersion estimation
    if meanvarmodel == None:
        alpha_dispersion=alpha_val
        alpha_dispersion_mat=alpha_val
    else:
        #alpha_dispersion_mat=np.matrix(sk.MAP_sgrna_dispersion_estimate)
        alpha_dispersion_mat=np.matrix(meanvarmodel.get_glm_dispersion(normalized_datak,returnalpha=True))
        #logging.info(normalized_datak)
        #logging.info(alpha_dispersion_mat)
        alpha_dispersion=alpha_dispersion_mat.getA1()
    stored_alpha_dispersion=copy.copy(alpha_dispersion)

    # logging.info(stored_alpha_dispersion)
    # initial guess of beta
    baseline_sample_matrix=np.log(sk.nb_count[basesampleid,:])-np.log(size_mat[basesampleid,:])
    if restart==True or len(sk.beta_estimate)==0:
        beta_vec1=np.mean(baseline_sample_matrix,axis=0)
        beta_vec2=np.log(sk.nb_count[1:,:])-np.log(size_mat[1:,:])-beta_vec1
        beta1_meanval=np.mean(beta_vec2,axis=1)
        beta1_es_mat=linalg.inv(design_mat.getT()*design_mat+alpha_val*np.matrix(np.identity(design_mat.shape[1])))*design_mat.getT()*(beta1_meanval)
        beta_init_mat=np.vstack((beta_vec1.getT(),beta1_es_mat))
        eff_list=[1]*n # a list of size nsgRNA
        sk.eff_estimate=eff_list
    else:
        beta_init_mat=np.matrix([[x] for x in sk.beta_estimate])
        eff_list=sk.eff_estimate

    if removeoutliers==True:
        eff_index=[i for i,j in enumerate(eff_list+[1]*nbeta1) if j==1]

        converting_matrix=np.identity(sk.nb_count.shape[0]*sk.nb_count.shape[1])
        efficieent_grna_index=[i for i,x in enumerate(eff_list*nallsample) if x==1]
        converting_matrix=converting_matrix[np.ix_(efficieent_grna_index,)]

        full_datak=copy.copy(datak)
        full_extdesign_mat=copy.copy(extdesign_mat)
        full_size_vec=copy.copy(size_vec)
        full_beta_init_mat=copy.copy(beta_init_mat)
        full_beta_prior_inverse_variance=copy.copy(beta_prior_inverse_variance)
        if PPI_prior==True:
            full_beta_prior_mean=copy.copy(beta_prior_mean)
            beta_prior_mean=beta_prior_mean[eff_index,:]

        beta_prior_inverse_variance=[beta_prior_inverse_variance[i] for i in eff_index]
        beta_init_mat=beta_init_mat[eff_index,:]
        datak=converting_matrix*datak
        extdesign_mat=(converting_matrix*extdesign_mat)[:,eff_index]
        size_vec=converting_matrix*size_vec

        alpha_dispersion=[alpha_dispersion[i] for i in range(len(alpha_dispersion)) if eff_list[i%len(eff_list)]==1]

    n_iter=1
    beta1_se_mat=None
    beta1_new=None
    while(True):
        if PPI_prior==False:
            logmu_estimate=extdesign_mat*beta_init_mat # (nsample*nsgrna)*1 matrix
            mu_estimate=(np.multiply(np.exp(logmu_estimate),size_vec))
            sgrna_residule=(datak-mu_estimate)
        else:
            # beta=beta_prior + beta_sec
            # mu=mu_prior * mu_sec
            # alpha_sec=alpha+(1/mu)-(1/mu_sec)
            # K~N_B(mu,alpha)
            # K_sec(==K/mu_prior)~N_B(mu_sec,alpha_sec)
            if n_iter==1:
                beta_init_mat=beta_init_mat-beta_prior_mean
            logmu_prior_estimate=extdesign_mat*beta_prior_mean
            mu_prior_estimate=np.exp(logmu_prior_estimate)

            logmu_estimate=extdesign_mat*beta_init_mat
            mu_estimate=(np.multiply(np.exp(logmu_estimate),size_vec))
            sgrna_residule=((datak/mu_prior_estimate)-mu_estimate)

            mu_prior_list=[i[0] for i in mu_prior_estimate.tolist()]
            mu_sec_list=[i[0] for i in mu_estimate.tolist()]
            mu_list=[mu_prior_list[i]*mu_sec_list[i] for i in range(len(mu_prior_list))]
            inverse_mu_sec_list=[(1/(i+10**(-10))) for i in mu_sec_list]
            inverse_mu_list=[1/(i+10**(-10)) for i in mu_list]
            alpha_sec_dispersion=[stored_alpha_dispersion[i]+inverse_mu_list[i]-inverse_mu_sec_list[i] for i in range(len(mu_list))]
            alpha_sec_dispersion_mat=(np.matrix(alpha_sec_dispersion)).reshape(len(alpha_sec_dispersion),1)
            alpha_dispersion=alpha_sec_dispersion

        z_estimate=sgrna_residule/mu_estimate + logmu_estimate

        w_list=mu_estimate.getA1()
        w_list_ele=1.0/((1.0/w_list)+alpha_dispersion)

        w_matrix=np.diag(w_list_ele)
        xwx_mat=extdesign_mat.getT()*w_matrix*extdesign_mat+np.diag(beta_prior_inverse_variance)
        xwx_inv=linalg.inv(xwx_mat)
        beta_new=xwx_inv*extdesign_mat.getT()*w_matrix*z_estimate

        beta_diff=beta_new-beta_init_mat
        if np.isnan(np.sum(beta_new)):
            break
        else:
            beta_init_mat=beta_new
        n_iter+=1
        diffval=np.sum(beta_diff.getA1()[n:]**2)
        absval=np.sum(beta_new.getA1()[n:]**2)
        if abs(absval)<1e-9:
            absval=1.0
        difffrac=diffval/absval

        if difffrac<diff_cutoff  or n_iter>n_max_init:
            break

    if PPI_prior==True:
        beta_init_mat=beta_init_mat+beta_prior_mean
        logmu_estimate=extdesign_mat*beta_init_mat
        mu_estimate=(np.multiply(np.exp(logmu_estimate),size_vec))
        sgrna_residule=(datak-mu_estimate)
        if removeoutliers==False:
            alpha_dispersion=stored_alpha_dispersion
        else:
            alpha_dispersion=[stored_alpha_dispersion[i] for i in range(len(stored_alpha_dispersion)) if eff_list[i%len(eff_list)]==1]
        w_list_ele=1.0/((1.0/mu_estimate.getA1())+alpha_dispersion)
        w_matrix=np.diag(w_list_ele)
        xwx_mat=extdesign_mat.getT()*w_matrix*extdesign_mat+np.diag(beta_prior_inverse_variance)
        xwx_inv=linalg.inv(xwx_mat)

    hat_matrix=np.sqrt(w_matrix)*extdesign_mat*xwx_inv*extdesign_mat.getT()*np.sqrt(w_matrix)
    v=extdesign_mat.shape[0]-(2*hat_matrix-hat_matrix*hat_matrix.getT()).trace()
    temp=((sgrna_residule/mu_estimate).getT()*(sgrna_residule/mu_estimate))/(v)
    beta1_se_mat=xwx_inv*extdesign_mat.getT()*w_matrix*w_matrix*extdesign_mat*xwx_inv
    beta_se_val=(temp.getA1()[0]**(0.5))*np.sqrt(np.diag(beta1_se_mat))

    if removeoutliers==True:
        beta_new_zscore_list=(beta_init_mat.getA1()/beta_se_val)[len([i for i in eff_list if i==1]):]
    else:
        beta_new_zscore_list=(beta_init_mat.getA1()/beta_se_val)[n:]
    beta_new_pval=norm.sf(np.abs(beta_new_zscore_list))*2

    if removeoutliers==True:
        eff_index=[[i,j] for i,j in enumerate(eff_index)]
        for k in eff_index:
            full_beta_init_mat[k[1]]=beta_init_mat[k[0]]
        sk.beta_estimate=full_beta_init_mat.getA1()
        logmu_estimate=full_extdesign_mat*full_beta_init_mat
        mu_estimate=(np.multiply(np.exp(logmu_estimate),full_size_vec))
    else:
        sk.beta_estimate=beta_init_mat.getA1()
        logmu_estimate=extdesign_mat*beta_init_mat
        mu_estimate=(np.multiply(np.exp(logmu_estimate),size_vec))


    #sk.beta_estimate=beta_init_mat.getA1()
    if removeoutliers==True:
        sk.beta_se_val=beta_se_val[len([i for i in eff_list if i==1]):]
    else:
        sk.beta_se_val=beta_se_val[n:]
    sk.beta_pval=beta_new_pval
    sk.beta_zscore=beta_new_zscore_list
    sk.beta_pval_pos=norm.sf(beta_new_zscore_list)
    sk.beta_pval_neg=norm.cdf(beta_new_zscore_list)

    sk.mu_estimate=mu_estimate
    sk.dispersion_estimate=stored_alpha_dispersion

    if meanvarmodel != None:
        sk.loglikelihood=getloglikelihood3(sk.sgrna_kvalue,sk.mu_estimate,stored_alpha_dispersion,sumup=False,log=True)
