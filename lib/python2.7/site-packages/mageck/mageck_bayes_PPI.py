from __future__ import print_function
import os
import sys
import math
import numpy as np
import pickle
from scipy import stats
import matplotlib.pyplot as pylab
import operator
from collections import Counter,defaultdict
import logging


def PPI_import_string_9(self):
    try:
        self.PPI_data=pickle.load(open("/Users/chen-haochen/Documents/Data_storage/CRISPR/pathway/string_9/string_v9_network.p",'rb'))
    except:
        self.PPI_data=pickle.load(open("/n/xiaoleliu_lab/chenhaochen/mageck_bayes/string_v9_network.p",'rb'))

def PPI_weighting_rewiring(self):
    weighting_record=[]
    for gene, neighboring_gene_weighting in self.PPI_data.items():
        weighting_record+=list(neighboring_gene_weighting.values())
    weighting_correction_ratio=np.sum(weighting_record)/len(weighting_record)

    for gene, neighboring_gene_weighting in self.PPI_data.items():
        for neighboring_gene, weighting in neighboring_gene_weighting.items():
            neighboring_gene_weighting[neighboring_gene]=weighting/weighting_correction_ratio

def PPI_coverge_diagnosis(self):
    PPI_genes=self.PPI_data.keys()
    crispr_genes=self.allgenedict.keys()
    PPI_genes=[i.upper() for i in PPI_genes]
    crispr_genes=[i.upper() for i in crispr_genes]
    logging.info("Number of PPI genes: %s" %str(len(PPI_genes)))
    logging.info("Number of CRISPR genes: %s" %str(len(crispr_genes)))
    logging.info("Number of overlapped genes: %s" %str(len([i for i in PPI_genes if i in crispr_genes])))
    if len([i for i in PPI_genes if i in crispr_genes])<3000:
        return False
    else:
        return True

def PPI_network_plotting(self,plotting=False,self_adjust=False,constant=0):
    beta1_estimate=defaultdict(list)
    beta1_se_mat=defaultdict(list)
    for (tgid,tginst) in self.allgenedict.items():
        tgid=tgid.upper()
        beta1_estimate[tgid]=tginst.beta_estimate.tolist()[tginst.nb_count.shape[1]:]
        beta1_se_mat[tgid]=tginst.beta_se_val
        #beta1_se_mat[tgid]=tginst.beta_se_val.tolist()

    true_central_beta=[]
    predicted_central_beta=[]
    gene_set=set(list(beta1_estimate.keys()))
    outliers_record=[]
    for gene,beta_value in beta1_estimate.items():
        #sgRNA_number=Counter(self.allgenedict[gene].eff_estimate)[1]
        if Counter([math.isnan(i) for i in beta_value])[True]==0:
            sgRNA_number=self.allgenedict[gene].nb_count.shape[1]
            original_surrounding_beta=[] # beta value of surrounding genes
            original_surrounding_beta_variance=[] # variance of beta of surrounding genes
            weighted_weighting=[] # the weighting given by the network data base
            for neighboring_gene, weighting in self.PPI_data[gene.upper()].items():
                if neighboring_gene in gene_set and Counter([math.isnan(i) for i in beta1_estimate[neighboring_gene]])[True]==0:
                    original_surrounding_beta.append(beta1_estimate[neighboring_gene])
                    original_surrounding_beta_variance.append([i**2 for i in beta1_se_mat[neighboring_gene]])
                    weighted_weighting.append(weighting**self.PPI_weighting)
            if len(weighted_weighting)>0:
                #logging.info(gene)
                #logging.info(beta_value)
                original_surrounding_beta=np.matrix(np.vstack(original_surrounding_beta))
                original_surrounding_beta_variance=np.matrix(np.vstack(original_surrounding_beta_variance))
                weighted_weighting=np.matrix(np.vstack(weighted_weighting)).T
                weighted_weighting_sum=np.sum(weighted_weighting,axis=1)
                weighted_average=((weighted_weighting*original_surrounding_beta)/(weighted_weighting_sum+constant))
                #weighted_average=((weighted_weighting*original_surrounding_beta))
                #logging.info(weighted_average.shape)

                #extended_weighted_average=np.tile(weighted_average,(original_surrounding_beta.shape[0],weighted_average.shape[1]))
                #weighted_residue=(original_surrounding_beta-extended_weighted_average)
                #weighted_residue_square=np.square(weighted_residue)
                #self.allgenedict[gene].prior_variance=((weighted_weighting*(weighted_residue_square+original_surrounding_beta_variance))/weighted_weighting_sum).tolist()[0]
                if self_adjust==True:
                    self.allgenedict[gene].prior_mean=weighted_average.tolist()[0]
                    self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance
                if Counter([math.isnan(i) for i in beta_value])[True]==0:
                    true_central_beta.append(beta_value)
                    predicted_central_beta.append(weighted_average.tolist()[0])
            else:
                outliers_record.append(gene)
                if self_adjust==True:
                    self.allgenedict[gene].prior_mean=[0]*len(beta_value)
                    self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance

        else:
            outliers_record.append(gene)
            if self_adjust==True:
                self.allgenedict[gene].prior_mean=[0]*len(beta_value)
                self.allgenedict[gene].prior_variance=self.non_PPI_beta_prior_variance

    true_central_beta=np.vstack(true_central_beta)
    predicted_central_beta=np.vstack(predicted_central_beta)

    if true_central_beta.shape!=predicted_central_beta.shape:
        sys.exit(1)

    beta_correlation_record=[]
    for k in range(true_central_beta.shape[1]):
        pearson_value=stats.pearsonr(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist())
        regression_value=np.polyfit(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist(),1)
        fit_fn = np.poly1d(regression_value)
        beta_correlation_record.append([pearson_value,regression_value])
        print(pearson_value)
        if plotting==True:
            #os.chdir(self.output_directory)
            pylab.scatter(predicted_central_beta[:,k].tolist(),true_central_beta[:,k].tolist(),color="r",s=0.1)
            pylab.plot(predicted_central_beta[:,k].tolist(), fit_fn(predicted_central_beta[:,k].tolist()), '--k')
            print ("Weighting ratio is {},y={}*x + {}".format(self.PPI_weighting,regression_value[0],regression_value[1]))
            pylab.xlabel("Weighted average of neighboring genes")
            pylab.ylabel("Center gene")
            pylab.savefig("{}_PPI_fitting_weighting_{}_beta1_{}_constant_{}.png".format(self.output_prefix,str(self.PPI_weighting),str(k),str(constant)))
            pylab.close()

    if self_adjust==False or plotting==False:
        return beta_correlation_record


def PPI_network_weighting_searching(self,constant=0):
    weighting_record=[]
    for weighting_ratio in [0+float(k)/5 for k in range(12)]:
        self.PPI_weighting=weighting_ratio
        beta_correlation_record=PPI_network_plotting(self,constant=constant)
        print(beta_correlation_record)
        weighting_record.append([weighting_ratio,beta_correlation_record[0][0][0]])

    os.chdir(self.output_directory)
    pylab.plot([i[0] for i in weighting_record],[i[1] for i in weighting_record])
    pylab.xlabel("Weighting ratio")
    pylab.ylabel("Correlation coefficient")
    pylab.savefig("{}_Weighting ratio selection.png".format(self.output_prefix))
    pylab.close()
    
    print (weighting_record)
    weighting_record.sort(key=operator.itemgetter(1),reverse=True)
    self.PPI_weighting=weighting_record[0][0]
    PPI_network_plotting(self,plotting=True,self_adjust=True,constant=constant)

def PPI_main(self,final_evaluation=False):
    constant=3
    PPI_import_string_9(self)
    if final_evaluation==True:
        PPI_weighting_rewiring(self)
        beta_correlation_record=PPI_network_plotting(self,plotting=True,self_adjust=False,constant=constant)
        return beta_correlation_record
    else:
        self.PPI_diagnosis=PPI_coverge_diagnosis(self)
        logging.info(self.PPI_diagnosis)
        if self.PPI_diagnosis==True:
            PPI_weighting_rewiring(self)
            if self.PPI_weighting==None:
                PPI_network_weighting_searching(self,constant=constant)
            else:
                #for constant in [0,1,2,3,4,5,6,10,20,50]:
                #    logging.info(constant)
                PPI_network_plotting(self,plotting=True,self_adjust=True,constant=constant)
