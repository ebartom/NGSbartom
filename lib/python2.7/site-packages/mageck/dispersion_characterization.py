import os
import copy
import numpy as np
import operator
from scipy.stats import nbinom
from scipy.special import factorial,gamma
import scipy
import scipy.stats
from collections import defaultdict
from mleclassdef import *
from mledesignmat import *
import logging
import decimal
#from sympy.matrices import *

def getloglikelihood2(k_list,mu_list,alpha,sumup=False,log=True):
    '''
    Get the log likelihood estimation of NB, using the current estimation of beta and alpha
    '''
    # solution 1
    mu_sq=np.multiply(mu_list,mu_list)
    var_vec=mu_list+np.multiply(alpha, mu_sq)
    nb_p=np.divide(mu_list,var_vec)
    nb_r=np.divide(mu_sq,var_vec-mu_list)

    if log:
        logp=nbinom.logpmf(k_list,nb_r,nb_p)
    else:
        logp=nbinom.pmf(k_list,nb_r,nb_p)
    if np.isnan(np.sum(logp)):
        logp=np.where(np.isnan(logp),0,logp)
    #print("hi",np.sum(logp))
    if sumup:
        #print(np.sum(logp))
        return np.sum(logp)
    else:
        #pass
        return logp

    # solution 2
    #temp=[]
    #print(k_list,mu_list)
    #for i in range(len(k_list)):
    #    temp_1=gamma(k_list[i]+(1/alpha))
    #    temp_2=factorial(k_list[i])
    #    temp_3=gamma(1/alpha)
    #    temp_4=((mu_list[i]*alpha)/(1+(mu_list[i]*alpha)))**(k_list[i])
    #    temp_5=(1/(1+mu_list[i]*alpha))**(1/alpha)
        #print(temp_1,temp_2,temp_3,temp_4,temp_5)
    #    temp.append(np.log((temp_1/(temp_2*temp_3))*temp_4*temp_5))
    #print(np.sum(temp))

    # solution 3
    #for i in range(len(k_list)):
    #    nb_r=1/alpha
    #    nb_p=(mu_list[i]*alpha)/(1+(mu_list[i]*alpha))

def Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha,extended_design_mat):
    log_likelihood_sum=getloglikelihood2(k_list,mu_list,alpha,sumup=True,log=True)

    w_matrix_ele=[1.0/(1.0/mu_list[i]+alpha) for i in range(len(mu_list))]
    w_matrix=np.diag(w_matrix_ele)
    xwx_mat=extended_design_mat.getT()*w_matrix*extended_design_mat
    Cox_Reid_adjustment=(0.5)*np.log(np.linalg.det(xwx_mat))

    Cox_Reid_adjusted_log_likelihood_sum=log_likelihood_sum-Cox_Reid_adjustment
    return(Cox_Reid_adjusted_log_likelihood_sum)

def first_derivative_of_log_likelihood(mu_list,k_list,alpha):
    temp0=0
    for i in range(len(mu_list)):
        temp_2=(1/(alpha**2))*(np.log(1+mu_list[i]*alpha))
        temp_3=(k_list[i]-mu_list[i])/(alpha*(1+mu_list[i]*alpha))
        temp_1=0
        for j in range(int(k_list[i])):
            temp_1-=(1/(alpha*(1+alpha*j)))
        temp0+=(temp_1+temp_2+temp_3)
    return(temp0)

def adjugate(matrix):
    C = np.zeros(matrix.shape)
    nrows, ncols = C.shape
    # Loop to calculate Cofactor
    for row in range(nrows):
        for col in range(ncols):
            minor = matrix[np.array(list(range(row))+list(range(row+1,nrows)))[:,np.newaxis],
                           np.array(list(range(col))+list(range(col+1,ncols)))]
            C[row, col] = (-1)**(row+col) * np.linalg.det(minor)
    return C.transpose()

def adjusted_profile(mu_list,alpha):
    temp_4=0.5*alpha
    #logging.info(mu_list)
    #logging.info(alpha)
    w_matrix_ele=[1.0/(1.0/mu_list[i]+alpha) for i in range(len(mu_list))]
    w_matrix=np.diag(w_matrix_ele)
    adjugate_w=adjugate(w_matrix)
    temp_4=temp_4*((adjugate_w*w_matrix*w_matrix).trace())/np.linalg.det(w_matrix)
    return(temp_4)

def sgrna_wide_dispersion_estimation_MAP(tginst):
    kappa=0.1
    gene_alpha=[]

    sgrna_kvalue=tginst.sgrna_kvalue.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
    mu_estimate=tginst.mu_estimate.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
    for i in range(len(tginst.sgrnaid)):
        alpha_record=[]
        mu_list=mu_estimate[i,].tolist()[0]
        k_list=sgrna_kvalue[i,].tolist()[0]

        #naive_alpha_record=[]
        #for k in range(len(mu_list)):
        #    print(mu_list[k])
        #    naive_alpha_record.append(((mu_list[k]-k_list[k])**2-(mu_list[k]))/(mu_list[k])**2)
        #naive_dispersion=np.mean(naive_alpha_record)

        naive_dispersion=0.2
        alpha_record.append(naive_dispersion)
        while len(alpha_record)<=2 or abs(alpha_record[-1]-alpha_record[-2])>0.005:
            #w_matrix_ele=[1.0/(1.0/mu_list[i]+alpha_record[-1]) for i in range(len(mu_list))]
            #w_matrix=np.diag(w_matrix_ele)
            direction=alpha_record[-1]*first_derivative_of_log_likelihood(mu_list,k_list,alpha_record[-1])+(adjusted_profile(mu_list,alpha_record[-1]))

            alpha_update=np.log(alpha_record[-1])+kappa*direction
            alpha_record.append(np.exp(alpha_update))
            if len(alpha_record)>20 or alpha_record[-1]<0:
                break
        #print(len(alpha_record),alpha_record)
        x=([i+k*tginst.nb_count.shape[1] for k in range(tginst.nb_count.shape[0])])
        grna_design_mat=tginst.extended_design_mat[x,:]
        grna_design_mat=grna_design_mat[:,[i,tginst.nb_count.shape[1]]]
        #for alpha_temp in [alpha_record[0],alpha_record[-1]]:
        #    print(alpha_temp)
        #for g in [0.01,0.1,1,2,4,8,16,32,64,108]:
        #    alpha_temp=g*alpha_record[-1]
        print(Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha_record[-1],grna_design_mat))
        print(Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha_record[-1]-0.01,grna_design_mat))
        print(Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha_record[-1]+0.01,grna_design_mat))
        #logging.info(len(alpha_record))
        gene_alpha.append(max(0.000000001,alpha_record[-1]))
    tginst.MAP_sgrna_dispersion_estimate=gene_alpha*(tginst.nb_count.shape[0])
    logging.info(gene_alpha)

def sgrna_wide_dispersion_estimation_MAP_v2(tginst,design_matrix):
    ratio=0.4
    gene_alpha=[]

    #Tracer()()
    if tginst.sgrna_kvalue.shape[0] == tginst.nb_count.shape[0]*tginst.nb_count.shape[1]:
      sgrna_kvalue=tginst.sgrna_kvalue.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
      mu_estimate=tginst.mu_estimate.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
    else:
      sgrna_kvalue=tginst.sgrna_kvalue[:tginst.nb_count.shape[0]*tginst.nb_count.shape[1],:].reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
      mu_estimate=tginst.mu_estimate[:tginst.nb_count.shape[0]*tginst.nb_count.shape[1],:].reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]).T
    for i in range(len(tginst.sgrnaid)):
        step_size=3
        alpha_record=[]
        mu_list=mu_estimate[i,].tolist()[0]
        k_list=sgrna_kvalue[i,].tolist()[0]
        mu_list=[min(10**10,i) for i in mu_list]

        mean_mu=np.mean(mu_list)
        variance=np.sum([(mu_list[i]-k_list[i])**2 for i in range(len(mu_list))])/(len(mu_list)-1)

        naive_dispersion=max(0.000000001,(variance-mean_mu)/(mean_mu**2))
        alpha_record=[0,naive_dispersion]
        grna_design_mat=design_matrix

        while abs(alpha_record[-1]-alpha_record[-2])>0.005:
            slope=alpha_record[-1]*first_derivative_of_log_likelihood(mu_list,k_list,alpha_record[-1])+(adjusted_profile(mu_list,alpha_record[-1]))

            if slope>0:
                direction=1
            else:
                direction=-1

            #temp_2=Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,w_matrix,alpha_record[-1],grna_design_mat)
            temp_2=Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha_record[-1],grna_design_mat)
            while True:
                alpha_temp=np.exp(np.log(alpha_record[-1])+step_size*direction)
                temp_1=Cox_Reid_adjusted_profile_likelihood(k_list,mu_list,alpha_temp,grna_design_mat)
                if temp_1>temp_2+ratio*slope*step_size or step_size<0.0001:
                    break
                else:
                    step_size=step_size*0.6

            alpha_record.append(max(0.000000001,alpha_temp))
            if len(alpha_record)>20:
                break

        gene_alpha.append(alpha_record[-1])

    #logging.info(gene_alpha)
    tginst.MAP_sgrna_dispersion_estimate=gene_alpha*(tginst.nb_count.shape[0])
    #logging.info(sgrna_kvalue)
    #logging.info(mu_estimate)
    #logging.info(tginst.MAP_sgrna_dispersion_estimate)
