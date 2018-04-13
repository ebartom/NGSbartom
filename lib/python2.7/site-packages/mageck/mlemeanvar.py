'''
Modeling mean and variance
'''

from __future__ import print_function

import re
import sys
import scipy
from scipy.stats import nbinom
from scipy.stats import linregress
from scipy.optimize import curve_fit
import random
import math
import numpy as np
import numpy.linalg as linalg
import logging


# debug
try:
  from IPython.core.debugger import Tracer
except:
  pass

from mageck.mleclassdef import *
from mageck.mlelowess import *


def func(x,a,b):
    return (a/x)+b

class MeanVarModel:
  '''
  The Mean and Variance model class
  '''
  # variables for storing k (read counts) and residule
  list_k=None
  list_res=None
  # for lowess regression
  list_lowess=None
  lowess_search={} # used to identify the lowess regression values of log(variance+0.01)
  # for linear regression
  lm_intercept=0.0
  lm_coeff=0.0
  # for generalized linear regression
  glm_a0=0.0
  glm_a1=0.0
  # member functions
  def get_mean_var_residule(self,allgenedict):
    '''
    Getting the read counts and regression residules for each sgRNA
    part of the modeling the mean and variance
    '''
    list_k=[]
    list_res=[]
    for (gid,gsk) in allgenedict.iteritems():
      nsg=len(gsk.nb_count[0])
      nsample=len(gsk.nb_count)
      if len(gsk.sgrna_kvalue)>0:
        sg_k=[x[0] for x in gsk.sgrna_kvalue.tolist()]
        sg_residule=[x[0] for x in gsk.sgrna_residule.tolist()]
        if len(sg_k)>=nsg*nsample:
          list_k+=sg_k[:(nsg*nsample)]
          list_res+=sg_residule[:(nsg*nsample)]
    np_k=np.array(list_k)
    np_res=np.array(list_res)
    self.list_k=np_k
    self.list_res=np_res
  #
  def save_k_residule_to_file(self,filename):
    # save the k and residule to file
    fileoh=open(filename,'w')
    for i in range(len(self.list_k)):
      print('\t'.join([str(self.list_k[i]),str(self.list_res[i])]),file=fileoh)
    fileoh.close()
  def model_mean_var_by_lowess(self):
    '''
    Modeling the mean and variance by LOWESS regression
    Too slow, need some changes in the future
    '''
    k_log=np.log(self.list_k)
    var_log=np.log(np.square(self.list_res)+0.01)
    self.list_lowess=lowess(k_log,var_log)
    # save to a dictionary
    for i in range(len(self.list_k)):
      self.lowess_search[self.list_k[i]]=self.list_lowess[i]
  #
  def get_lowess_var(self,klist,returnalpha=False):
    '''
    Return the fitted values of variance.
    If returnalpha=True, return the alpha value (var=mean+alpha*mean^2)
    '''
    if self.list_lowess==None:
      raise ValueError('Lowess regression has not been done yet.')
    for k in klist:
      if k not in self.lowess_search:
        raise ValueError('the key value is not in the dictionary.')
    varvalue=np.array([self.lowess_search[k] for k in klist])
    varvalue=np.exp(varvalue)-0.01
    varvalue=np.array([ (lambda x: x if x>=0.01 else 0.01)(t) for t in varvalue])
    # 
    if returnalpha:
      #alphavalue=(varvalue-k)/(k**2)
      alphavalue=[(varvalue[i]-klist[i])/(klist[i]**2) for i in range(len(klist))]
      alphavalue=np.array([ (lambda x: x if x>=0.01 else 0.01)(t) for t in alphavalue])
      return alphavalue
    else:
      return varvalue
  def model_mean_var_by_lm(self):
    '''
    Modeling the mean and variance by linear regression
    '''
    k_log=np.log(self.list_k)
    var_log=np.log(np.square(self.list_res)+0.01)
    # remove those with too low variance
    k_log2=np.array([k_log[i] for i in range(len(var_log)) if var_log[i]>(-1)])
    var_log2=np.array([var_log[i] for i in range(len(var_log)) if var_log[i]>(-1)])
    if len(k_log2)>20:
      (slope,intercept,r_value,p_value,std_err)=linregress(k_log2,var_log2)
    else:
      (slope,intercept,r_value,p_value,std_err)=linregress(k_log,var_log)
    self.lm_intercept=intercept
    self.lm_coeff=slope
    # Tracer()()
    logging.info('Linear regression: y='+str(slope)+'x+'+str(intercept))
    if np.isnan(slope) or np.isnan(intercept):
      logging.error('Nan values for linear regression')
  
  def get_lm_var(self,klist,returnalpha=False):
    '''
    Return the fitted values of variance.
    If returnalpha=True, return the alpha value (var=mean+alpha*mean^2)
    '''
    kls=(klist)
    k_log=np.log(kls)
    varvalue=k_log*self.lm_coeff+self.lm_intercept
    varvalue=np.exp(varvalue)-0.01
    #varvalue=np.array([ (lambda x: x if x>=0.01 else 0.01)(t) for t in varvalue])
    varvalue=np.where(varvalue>0.01,varvalue,0.01)
    # set up the lower bound
    th_count=10.0
    var_t=np.log(th_count)*self.lm_coeff+self.lm_intercept
    var_t=np.exp(var_t)-0.01
    if var_t<0.01:
      var_t=0.01
    alpha_t=(var_t-th_count)/(th_count*th_count)
    if alpha_t<1e-2:
      alpha_t=1e-2
    # 
    if returnalpha:
      #alphavalue=[(varvalue[i]-kls[i])/(kls[i]**2) for i in range(len(kls))]
      #alphavalue=np.array([ (lambda x: x if x>=0.01 else 0.01)(t) for t in alphavalue])
      alphavalue=np.divide((varvalue-kls),np.multiply(kls,kls))
      alphavalue2=np.where(alphavalue>alpha_t,alphavalue,alpha_t)
      return alphavalue2
    else:
      return varvalue
   

  def model_mean_disp_by_glm(self,allgenedict,output_prefix,size_f):
      '''
      Modeling the mean and dispersion by linear regression
      '''
      list_k=[]
      list_dispersion=[]
      '''
      for (gid,gsk) in allgenedict.iteritems():
          nsg=len(gsk.nb_count[0])
          nsample=len(gsk.nb_count)
          if len(gsk.sgrna_kvalue)>0:
              if gsk.MAP_sgrna_dispersion_estimate!=None:
                  sg_k=[x[0] for x in gsk.sgrna_kvalue.tolist()]
                  sg_dispersion=gsk.MAP_sgrna_dispersion_estimate
                  if len(sg_k)>=nsg*nsample:
                      list_k+=sg_k[:(nsg*nsample)]
                      list_dispersion+=sg_dispersion[:(nsg*nsample)]
      '''
      inverse_size_f=[1/i for i in size_f]
      for (gid,gsk) in allgenedict.items():
          nallsample=gsk.nb_count.shape[0]
          n=gsk.nb_count.shape[1]
          if gsk.MAP_sgrna_dispersion_estimate!=None:
              sg_dispersion=gsk.MAP_sgrna_dispersion_estimate
              sg_k=[x[0] for x in gsk.sgrna_kvalue.tolist()]
              list_dispersion+=sg_dispersion[:len(gsk.sgrnaid)]
              for i in range(len(gsk.sgrnaid)):
                  normalized_k_mean=np.mean(np.multiply(inverse_size_f,[sg_k[i+k*n] for k in range(nallsample)]))
                  #normalized_k_mean=np.mean(np.multiply(inverse_size_f,sg_k[i*nallsample:(i+1)*nallsample]))
                  list_k.append(normalized_k_mean)


      #combined_list=[[list_k[i],list_dispersion[i]] for i in range(len(list_k)) if math.isnan(list_k[i])==False and math.isnan(list_dispersion[i])==False]
      #combined_list=[[np.log(list_k[i]),np.log(list_dispersion[i])] for i in range(len(list_k)) if math.isnan(list_k[i])==False and math.isnan(list_dispersion[i])==False]
      combined_list=[[list_k[i],list_dispersion[i]] for i in range(len(list_k)) if math.isnan(list_k[i])==False and list_dispersion[i]<10 and math.isnan(list_dispersion[i])==False]
      xdata=np.array([i[0] for i in combined_list])
      ydata=np.array([i[1] for i in combined_list])
      #Tracer()()
      popt, pcov = curve_fit(func, xdata, ydata)
      logging.info(popt)
      logging.info(pcov)
      self.glm_a0=popt[0]
      self.glm_a1=popt[1]
      # pylab
      '''
      xdemo = np.linspace(0, 2000, 500)
      ydemo = func(xdemo,popt[0],popt[1])
      pylab.scatter([i[0] for i in combined_list],[i[1] for i in combined_list],s=0.01)
      pylab.plot(xdemo,ydemo,'-')
      pylab.savefig("mean_dispersion_%s" %output_prefix)
      pylab.close()
      '''
      # 
      '''
      matrix_list_dispersion=np.matrix(ydata).reshape(len(ydata),1)
      matrix_inverse_k_list=np.matrix([[i,1] for i in xdata]).reshape(len(ydata),2)
      gamma_model = smf.GLM(matrix_list_dispersion,matrix_inverse_k_list,family=sm.families.Gamma()).fit()
      print(gamma_model.summary())
      para=gamma_model.params
      self.glm_a0=para[1]
      self.glm_a1=para[0]
      '''
  def get_glm_dispersion(self,klist,returnalpha=False):
      kls=np.array(klist)
      #logging.info(kls)
      #dispvalue=(self.glm_a1*kls+self.glm_a0)
      #dispvalue=np.exp(self.glm_a1*np.log(kls)+self.glm_a0)
      #logging.info(self.glm_a0)
      #logging.info(self.glm_a1)
      dispvalue=func(kls,self.glm_a0,self.glm_a1)
      return(dispvalue)



