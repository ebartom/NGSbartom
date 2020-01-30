import os
import pickle
import numpy as np
import matplotlib.pyplot as pylab
import operator
from multiprocessing import Pool, Manager
import glob
from collections import Counter,defaultdict
from scipy.stats import norm
import logging

def bayes_selection_constnat_optimization(self):
    removal_ratio=0.9

    likelihood_list=[]
    modified_likelihood_list=[]
    sgRNA_number=[] # total sgRNA_number

    for (tgid,tginst) in self.allgenedict.items():
        #probability=tginst.sgrna_probability
        likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
        likelihood_mean= np.matrix(np.mean(likelihood,axis=0))
        #modified_likelihood_list+=[likelihood_mean.tolist()[0][i]-probability[i] for i in range(len(probability))]
        modified_likelihood_list+=likelihood_mean.tolist()[0]
        #probability_list+=probability
        likelihood_list+=(likelihood_mean.tolist()[0])
        sgRNA_number.append(len(tginst.sgrnaid))

    # selection_constant iteration start point
    sgRNA_number=[[sgRNA_number,count] for sgRNA_number,count in Counter(sgRNA_number).items()]
    sgRNA_number.sort(key=operator.itemgetter(1),reverse=True)
    common_sgRNA_number=sgRNA_number[0][0]
    modified_likelihood_list.sort(reverse=True)
    selection_constant=int(-modified_likelihood_list[int(len(modified_likelihood_list)*removal_ratio)]/(np.log(int(common_sgRNA_number*removal_ratio))-np.log(int(common_sgRNA_number*removal_ratio)-1)))

    # log_precaculation
    #log_list=defaultdict(list)
    for i in range(1,max([i[0] for i in sgRNA_number])+1):
        self.log_list[i]=[np.log(k+1) for k in range(i)]

    # selection constant searching
    constant_percentage={} # to record the percenge of on-target gRNAs given selection constant
    while selection_constant not in list(constant_percentage.keys()):
        total_sgRNA=0
        total_on_target_sgRNA=0
        for (tgid,tginst) in self.allgenedict.items():
            likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
            likelihood_mean= np.matrix(np.mean(likelihood,axis=0)).tolist()
            temp=likelihood_mean[0]
            #probability=tginst.sgrna_probability
            #temp=[temp[i]-probability[i] for i in range(len(temp))]
            temp.sort(reverse=True)
            temp_accumulate=[np.sum([k for k in temp[:(i+1)]]) for i in range(len(temp))]
            modified_penalty=[selection_constant*i for i in self.log_list[len(tginst.sgrnaid)]]
            temp_accumulate_log_penalty=[sum(x) for x in zip(temp_accumulate,modified_penalty)]
            max_index=[i for i, j in enumerate(temp_accumulate_log_penalty) if j == max(temp_accumulate_log_penalty)]
            total_sgRNA+=len(tginst.sgrnaid)
            total_on_target_sgRNA+=(max_index[0]+1)
            if abs(total_sgRNA-3000)<10:
                break
        constant_percentage[selection_constant]=float(total_on_target_sgRNA)/total_sgRNA
        if float(total_on_target_sgRNA)/total_sgRNA>removal_ratio:
            selection_constant-=1
        else:
            selection_constant+=1
        if [i>removal_ratio for i in list(constant_percentage.values())]==[True]*len(list(constant_percentage.values())):
            pass
        elif [i>removal_ratio for i in list(constant_percentage.values())]==[False]*len(list(constant_percentage.values())):
            pass
        else:
            break

    selection_percentage=[[constant,abs(percentage-removal_ratio)] for constant,percentage in constant_percentage.items()]
    selection_percentage.sort(key=operator.itemgetter(1))
    selection_constant=selection_percentage[0][0]
    self.selection_constant=selection_constant
    logging.info(selection_constant)

def bayes_selection(tginst,log_list,selection_constant):
    #likelihood_mean=np.mean(tginst.loglikelihood[i])
    likelihood= np.matrix(tginst.loglikelihood.reshape(tginst.nb_count.shape[0],tginst.nb_count.shape[1]))
    likelihood_mean= np.matrix(np.mean(likelihood,axis=0)).tolist()[0]
    #temp=[[tginst.sgrnaid[i],likelihood_mean[0][i]-probability[i]] for i in range(len(tginst.sgrnaid))]
    #logging.info(tginst.loglikelihood)
    #logging.info(likelihood_mean)
    #for i in range(len(tginst.sgrnaid)):
    #    logging.info(tginst.sgrnaid[i])
    #    logging.info(tginst.loglikelihood[i])
    temp=[[tginst.sgrnaid[i],likelihood_mean[i]] for i in range(len(tginst.sgrnaid))]
    temp.sort(key=operator.itemgetter(1),reverse=True)
    #logging.info(temp)

    temp_accumulate=[np.sum([k[1] for k in temp[:(i+1)]]) for i in range(len(temp))]

    modified_penalty=[selection_constant*i for i in log_list[len(tginst.sgrnaid)]]
    temp_accumulate_log_penalty=[sum(x) for x in zip(temp_accumulate,modified_penalty)]
    #logging.info(temp_accumulate_log_penalty)
    #if max(temp_accumulate_log_penalty)<0:
    max_index=[i for i, j in enumerate(temp_accumulate_log_penalty) if j == max(temp_accumulate_log_penalty)]
    #else:
    #    max_index=[i for i, j in enumerate(temp_accumulate_log_penalty) if j>*max(temp_accumulate_log_penalty)]
    on_target_sgRNA=[temp[i][0] for i in range(max_index[-1]+1)]
    tginst.eff_estimate=[1]*len(tginst.sgrnaid)
    outliers_index=[orders for orders,sgRNA in enumerate(tginst.sgrnaid) if sgRNA not in on_target_sgRNA]
    for i in outliers_index:
        tginst.eff_estimate[i]=0

def background(self):
    sgRNA_beta=dict()
    sgRNA_probability=dict()
    os.chdir(self.output_directory)
    #os.chdir("/Users/chen-haochen/Documents/Data_storage/CRISPR/MAGECK/Tim/singlized_Tim.txt_mageck_bayes")
    data=pickle.load(open("self_pre_dispersion_fitting_round.p",'rb'))
    for (tgid,tginst) in data.allgenedict.items():
        sgRNA_beta[tgid]=tginst.beta_estimate[1]
    return(sgRNA_beta)
