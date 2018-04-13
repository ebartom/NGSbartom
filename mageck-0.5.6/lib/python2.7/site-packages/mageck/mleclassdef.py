'''
Class definition
'''

from __future__ import print_function
# from IPython.core.debugger import Tracer

class SimCaseSimple:
  prefix='sample1'
  # the beta parameters; beta0 is the base value (a list of double, size nsgRNA)
  beta0=[]
  # beta1 (a nsample*r) is the beta values of different conditions
  beta1=[0]
  # the efficient parameters; binary, 0 means it's not efficient and 1 means it's efficient
  isefficient=[]
  # NB parameters used to generate read counts
  #mu0=[]
  #mu1=[]
  mu_estimate=[]
  #var0=[]
  #var1=[]
  var=[]
  #nb_p0=[]
  #nb_p1=[]
  nb_p=[]
  #nb_r0=[]
  #nb_r1=[]
  nb_r=[]
  # actual read counts
  #nb_count0=[]
  #nb_count1=[]
  nb_count=[]
  # design matrix
  design_mat=[]
  # extended design matrix
  extended_design_mat=[]
  extended_design_mat_residule=[]
  #
  # estimated values
  beta_estimate=[]
  w_estimate=[]
  # p values
  beta_zscore=[]
  beta_pval=[]  # two-sided p values
  beta_pval_fdr=[]  # two-sided p values
  beta_pval_neg=[]  # one-sided p values
  beta_pval_neg_fdr=[]  # one-sided p values
  beta_pval_pos=[]  # one-sided p values
  beta_pval_pos_fdr=[]  # one-sided p values
  beta_permute_pval=[]
  beta_permute_pval_fdr=[]
  beta_permute_pval_neg=[]
  beta_permute_pval_neg_fdr=[]
  beta_permute_pval_pos=[]
  beta_permute_pval_pos_fdr=[]
  # sgRNA ids
  sgrnaid=[]
  # residules, used for calculating the mean-variance modeling
  sgrna_kvalue=[] # count matrix; is usually (2*(nsample)+1)*1
  sgrna_residule=[] # fitted residule; the same size as sgrna_kvalue
  # dispersion estimate (i.e., r, or alpha)
  dispersion_estimate=None
  # MAP estimate
  MAP_sgrna_dispersion_estimate=None
  MAP_gene_dispersion_estimate=None
  non_PPI_beta_prior_variance=None
  prior_mean=None
  prior_variance=None
  sgrna_probability=None
  loglikelihood=None
  def gene_to_printfield(self,onesided=False):
    '''
    Convert values to print field
    '''
    nsg=0
    if len(self.nb_count)>0:
      nsg=(self.nb_count.shape[1])
    # to field
    ret=[self.prefix,str(nsg)]
    for i in range(nsg,len(self.beta_estimate)):
      ret+=['{0:.5g}'.format(self.beta_estimate[i])]
      ret+=['{0:.5g}'.format(self.beta_zscore[i-nsg])]
      if onesided:
        # one sided test
        ret+=['{0:.5g}'.format(self.beta_pval_neg[i-nsg]), '{0:.5g}'.format(self.beta_pval_neg_fdr[i-nsg])]
        ret+=['{0:.5g}'.format(self.beta_pval_pos[i-nsg]), '{0:.5g}'.format(self.beta_pval_pos_fdr[i-nsg])]
        if len(self.beta_permute_pval)>0:
          ret+=['{0:.5g}'.format(self.beta_permute_pval_neg[i-nsg]), '{0:.5g}'.format(self.beta_permute_pval_neg_fdr[i-nsg])]
          ret+=['{0:.5g}'.format(self.beta_permute_pval_pos[i-nsg]), '{0:.5g}'.format(self.beta_permute_pval_pos_fdr[i-nsg])]
        else:
          ret+=['1.0','1.0','1.0','1.0']
      else:
        # two-sided test
        if len(self.beta_permute_pval)>0:
          ret+=['{0:.5g}'.format(self.beta_permute_pval[i-nsg]),'{0:.5g}'.format(self.beta_permute_pval_fdr[i-nsg])]
        else:
          ret+=['1.0','1.0']
        ret+=['{0:.5g}'.format(self.beta_pval[i-nsg]), '{0:.5g}'.format(self.beta_pval_fdr[i-nsg])]
    return ret


def decformat(x):
  return '{0:.3g}'.format(x)


def gene_fdr_correction(allgenedict, method):
  '''
  perform p value correction
  '''
  pvaluelabel_list=['beta_pval','beta_pval_neg','beta_pval_pos','beta_permute_pval','beta_permute_pval_neg','beta_permute_pval_pos']
  fdr_label_list=[x+'_fdr' for x in pvaluelabel_list]
  for ii in range(len(pvaluelabel_list)):
    pvaluemat_list=[]
    var_fdr_list=[]
    #whichp='beta_pval'
    #writep='beta_pval_fdr'
    whichp=pvaluelabel_list[ii]
    writep=fdr_label_list[ii]
    for (gene,ginst) in allgenedict.iteritems():
      tlist=getattr(ginst,whichp)
      pvaluemat_list+=[tlist]
    # 
    import numpy as np
    from mageck.fdr_calculation import pFDR
    pvaluemat=np.matrix(pvaluemat_list)
    pvaluemat_t=pvaluemat.getT()
    # row by row
    for cid in range(pvaluemat_t.shape[0]):
      p_vec=pvaluemat_t[cid,:].getA1().tolist()
      fdr_vec=pFDR(p_vec,method)
      #Tracer()()
      pvaluemat_t[cid,:]=np.array(fdr_vec)
    # set up the attribute
    gid=0
    for (gene,ginst) in allgenedict.iteritems():
      targetp=pvaluemat_t[:,gid].getA1().tolist()
      setattr(ginst,writep,targetp)
      gid+=1
  

