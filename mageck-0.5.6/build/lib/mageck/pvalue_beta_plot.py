import os
import scipy.stats as stats
import matplotlib.pyplot as pylab
import numpy as np
from scipy.stats import norm
import statsmodels.api as sm

def plot(file_name,negative_control_gRNAs=None,wald_only=False):
    data=open(file_name,'rb')
    short_file_name=file_name[:file_name.index(".gene_summary.txt")]
    data.readline()
    permute_p_value_list=[]
    wald_p_value_list=[]
    beta_value_list=[]

    if negative_control_gRNAs!=None:
        negative_control_permute_p_value_list=[]
        negative_control_wald_p_value_list=[]
        negative_control_beta_value_list=[]


    for line in data:
        elements=line.decode().strip().split("\t")
        if negative_control_gRNAs!=None and elements[0] in negative_control_gRNAs:
            negative_control_beta_value_list.append(float(elements[2]))
            if wald_only==True:
                negative_control_wald_p_value_list.append(float(elements[4]))
            else:
                negative_control_permute_p_value_list.append(float(elements[4]))
                negative_control_wald_p_value_list.append(float(elements[6]))
        else:
            beta_value_list.append(float(elements[2]))
            if wald_only==True:
                wald_p_value_list.append(float(elements[4]))
            else:
                permute_p_value_list.append(float(elements[4]))
                wald_p_value_list.append(float(elements[6]))
    beta_value_list=[x for x in beta_value_list if str(x) != 'nan' and abs(x)<3]
    wald_p_value_list=[x for x in wald_p_value_list if str(x) != 'nan']
    if negative_control_gRNAs!=None:
        negative_control_beta_value_list=[x for x in beta_value_list if str(x) != 'nan' and abs(x)<3]
        negative_control_wald_p_value_list=[x for x in wald_p_value_list if str(x) != 'nan']

    if wald_only!=True:
        permute_p_value_list=[x for x in permute_p_value_list if str(x) != 'nan']
        stats.probplot(permute_p_value_list, dist="uniform",plot=pylab)
        pylab.savefig("QQplot of permute_p value %s.png" %short_file_name)
        pylab.close()

    pylab.hist(beta_value_list,bins=1000)
    pylab.savefig("Hist of beta value %s.png" %short_file_name)
    pylab.close()

    #stats.probplot(wald_p_value_list, dist="uniform",plot=pylab)
    fig=sm.qqplot(np.array(wald_p_value_list),stats.uniform,fit=True, line='45')
    pylab.xlim(0,1)
    pylab.ylim(0,1)
    #fig.set_xlim(0,1)
    pylab.savefig("QQplot of wald_p value %s.png" %short_file_name)
    pylab.close()
    '''
    if negative_control_gRNAs!=None:
        pylab.hist(negative_control_beta_value_list,bins=1000)
        pylab.savefig("Hist of negative control beta value %s.png" %short_file_name)
        pylab.close()

        stats.probplot(negative_control_wald_p_value_list, dist="uniform",plot=pylab)
        pylab.savefig("QQplot of negative control wald_p value %s.png" %short_file_name)
        pylab.close()
    '''
