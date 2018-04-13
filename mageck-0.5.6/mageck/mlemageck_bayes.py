#!/usr/bin/env python
'''
MAGeCK MLE main entry
'''

from __future__ import print_function

import re
import sys
import random
import math
import logging
import numpy as np
from scipy.stats import norm
import scipy
from scipy.stats import nbinom
import numpy as np
import numpy.linalg as linalg
from multiprocessing import Pool, Manager
from collections import defaultdict
import pickle
import glob

from mleinstanceio import *
from mleem_bayes import iteratenbem
from mlemeanvar import MeanVarModel
from mageckCount import normalizeCounts
from bayes_selection import *
from dispersion_characterization import *
from mleargparse import *
from mageck_bayes_PPI import *
from pvalue_beta_plot import *
from mleclassdef import *
# from mageck_bayes_QC import *

def sgRNA_production(file):
    data=open(file,'rb')
    output=open("singlized_%s" %(file),'wb')

    output.write(data.readline())
    for line in data:
        elements=line.decode().strip().split("\t")
        insert=[elements[i] for i in ([0,0]+list(range(2,len(elements))))]
        output.write("%s\n" %"\t".join(insert))

def beta_non_PPI_prior_calculation(allgenedict):
    '''
    For zero-centered beta prior, this function calculates the variance
    '''
    temp_beta=[(sk.nb_count.shape[1],sk.beta_estimate) for gene,sk in allgenedict.items()]
    temp_beta1=[i[1][i[0]:] for i in temp_beta]
    array_beta1=np.asarray(temp_beta1)
    beta1_prior_variance=[]
    for column in range(array_beta1.shape[1]):
        beta1=array_beta1[:,column]
        normalized_beta1=[(i-np.mean(beta1)) for i in beta1]
        abs_beta1=[abs(k) for k in normalized_beta1]
        abs_beta1.sort()
        var_beta1=float(np.percentile(abs_beta1,95))/norm.ppf(0.975)
        beta1_prior_variance.append((var_beta1*2)**2)
    return(beta1_prior_variance)

def dispersion_reestimate_sup(input_variable):
    self,tgid_list=input_variable
    if len(tgid_list)>0:
        for (tgid,tginst) in self.allgenedict.items():
            if tgid in tgid_list:
                sgrna_wide_dispersion_estimation_MAP(tginst)

        temp={tgid:tginst for tgid,tginst in self.allgenedict.items() if tgid in tgid_list}
        pickle.dump(temp,file=open("%s_dispersion_reestimate.p" %(tgid_list[0]),'wb'))

def crisprseq_parseargs():
    """
    Parsing mageck arguments.
    """
    parser=argparse.ArgumentParser(description='mageck: performs sgRNA, gene and pathway analysis on CRISPR-Cas9 screening data.')
    # definition of sub commands
    subparser=parser.add_subparsers(help='commands to run mageck',dest='subcmd')
    parser.add_argument('-v', '--version',action='version',version='%(prog)s 0.5.1')

    subm_mle=subparser.add_parser('bayes',help='Perform bayes estimation of gene essentiality.')
    subm_mle.add_argument('-n','--output-prefix',default=None,help='The prefix of the output file(s).')
    subm_mle.add_argument('--genes-varmodeling',default="all",help='The number of genes for mean-variance modeling. Default all.')
    subm_mle.add_argument('--permutation-round',type=int,default=3,help='The rounds for permutation (interger). The permutation time is (# genes)*x for x rounds of permutation. Suggested value: 100 (may take longer time). Default 2.')
    subm_mle.add_argument('-i', '--include-samples', help='Specify the sample labels if the design matrix is not given by file in the --design-matrix option. Sample labels are separated by ",", and must match the labels in the count table.')
    subm_mle.add_argument('-b', '--beta-labels', help='Specify the labels of the variables (i.e., beta), if the design matrix is not given by file in the --design-matrix option. Should be separated by ",", and the number of labels must equal to (# columns of design matrix), including baseline labels. Default value: "bata_0,beta_1,beta_2,...".')
    subm_mle.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr), holm\'s method (holm), or pounds\'s method (pounds).')
    subm_mle.add_argument("-o","--outliers_removal",action='store_true',help="Speicify whehter you want to remove outliers and recalculate..")

    subm_mle.add_argument("--norm-method",choices=['none','median','total','control'],default='median',help='Method for normalization, including "none" (nonormalization), "median" (median normalization, default), "total" (normalization by total read counts), "control" (normalization by control sgRNAs specified by the --control-sgrna option).')
    #subm_mle.add_argument("--control-sgrna",help="A file contains control sgRNAs")

    subm_mle.add_argument("-p","--PPI_prior",action='store_true',help="Specify whether you want to incorporate PPI as prior")
    subm_mle.add_argument("-w","--PPI_weighting",type=float,help="The weighting used to calculate PPI prior. If not provided, iterations will be used.",default=None)
    subm_mle.add_argument("-e","--negative_control",help="The gene name of negative controls. The corresponding sgRNA will be viewed independently.",default=None)

    # required parameters
    subm_mle.add_argument('-k','--count_table',required=True,help='Provide a tab-separated count table. Each line in the table should include sgRNA name (1st column), target gene (2nd column) and read counts in each sample.')
    subm_mle.add_argument('-d','--design_matrix',required=True,help='Provide a design matrix, either a file name or a quoted string of the design matrix. For example, "1,1;1,0". The row of the design matrix must match the order of the samples in the count table (if --include-samples is not specified), or the order of the samples by the --include-samples option.')

    args=parser.parse_args()

    if args.subcmd == None:
        parser.print_help()
        sys.exit(0)
    if args.output_prefix==None:
        args.output_prefix="{}_mageck_bayes".format(args.count_table)
    if args.negative_control!=None:
        args.negative_control=args.negative_control.split(',')
        args.negative_control=[i.upper() for i in args.negative_control]

    return args

def postargs(args):
    '''
    post-processing of argument parsing
    '''
    # configure logging information
    logging.basicConfig(level=10,
        format='%(levelname)-5s @ %(asctime)s.%(msecs)03d: %(message)s ',
        datefmt='%a, %d %b %Y %H:%M:%S',
        # stream=sys.stderr,
        filename=args.output_prefix+'.log',
        filemode='w'
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s.%(msecs)03d: %(message)s ','%a, %d %b %Y %H:%M:%S')
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    logging.info('Parameters: '+' '.join(sys.argv))

    from mledesignmat import parse_designmat

    try:
        import scipy
        from scipy.stats import nbinom
    except ImportError:
        logging.error('Cannot find scipy (required for mle approach). Please check your scipy installation.')
        sys.exit(-1)
    try:
        import numpy as np
        import numpy.linalg as linalg
    except ImportError:
        logging.error('Cannot find numpy (required for mle approach). Please check your numpy installation.')
        sys.exit(-1)
    # parsing design matrix
    (desmat,sampleid,betalabel)=parse_designmat(args.design_matrix)
    args.design_matrix=desmat

    # parsing sample label
    if sampleid ==None:
        # design matrix is provided as a string
        if args.include_samples !=None:
            args.include_samples=args.include_samples.split(',')
            if len(args.include_samples) != desmat.shape[0]:
                logging.error('The number of samples in the --include-samples option do not match rows in design matrix.')
                sys.exit(-1)
        if args.beta_labels!=None:
            args.beta_labels=args.beta_labels.split(',')
            if len(args.beta_labels) != desmat.shape[1]:
                logging.error('The number of labels in the --beta-labels option do not match columns in design matrix.')
                sys.exit(-1)
    else:
        # design matrix is provided as file
        if args.include_samples !=None:
            logging.error('Sample labels are included in the design matrix file '+args.design_matrix+'. The --include-samples option should not be used.')
            sys.exit(0)
        if args.beta_labels!=None:
            logging.error('Beta labels are included in the design matrix file '+args.design_matrix+'. The --beta-labels option should not be used.')
            sys.exit(0)
        args.include_samples=sampleid
        args.beta_labels=betalabel
        if len(args.include_samples) != desmat.shape[0]:
            logging.error('The number of samples in the --include-samples option do not match rows in design matrix.')
            sys.exit(-1)
        if len(args.beta_labels) != desmat.shape[1]:
            logging.error('The number of labels in the --beta-labels option do not match columns in design matrix.')
            sys.exit(-1)
    # log design matrix and column, row labels
    logging.info('Design matrix:')
    for desmat_1line in str(desmat).split('\n'):
        logging.info(desmat_1line)
    if args.beta_labels != None:
        logging.info('Beta labels:'+','.join(args.beta_labels))
    if args.include_samples != None:
        logging.info('Included samples:'+','.join(args.include_samples))

    return args

class Mageck_Bayes():
    def __init__(self, options):
        self.adjust_method=options.adjust_method
        self.beta_labels=options.beta_labels
        self.count_table=options.count_table
        self.design_matrix=options.design_matrix
        self.genes_varmodeling=options.genes_varmodeling
        self.include_samples=options.include_samples
        self.output_prefix=options.output_prefix
        self.file_directory=os.getcwd()
        self.output_directory="%s/%s" %(os.getcwd(),self.output_prefix)
        if os.path.exists(self.output_directory)==False:
            os.makedirs(self.output_directory)
        self.permutation_round=options.permutation_round
        self.outliers_removal=options.outliers_removal
        self.PPI_prior=options.PPI_prior
        self.PPI_weighting=options.PPI_weighting
        self.PPI_diagnosis=None
        self.non_PPI_beta_prior_variance=None
        self.allgenedict=None
        self.size_f=None
        self.mrm=None
        self.log_list=defaultdict(list)
        self.selection_constant=None
        self.invalid_gRNA_dict=None
        self.negative_control=options.negative_control
        self.negative_control_gRNAs=None
        self.norm_method=options.norm_method

    def bayes_init(self):
        maxgene=np.inf
        if self.negative_control==None:
            self.allgenedict,self.invalid_gRNA_dict=read_gene_from_file(self.count_table,includesamples=self.include_samples,negative_control=self.negative_control)
        else:
            self.allgenedict,self.invalid_gRNA_dict,self.negative_control_gRNAs=read_gene_from_file(self.count_table,includesamples=self.include_samples,negative_control=self.negative_control)

        # calculate the size factor
        cttab_sel={}
        for (geneid,gk) in self.allgenedict.iteritems():
            sgid=gk.sgrnaid
            sgreadmat=gk.nb_count.getT().tolist()
            for i in range(len(sgid)):
                cttab_sel[sgid[i]]=sgreadmat[i]
        if hasattr(self,'norm_method'):
            if self.norm_method!='none':
                self.size_f=normalizeCounts(cttab_sel,method=self.norm_method,returnfactor=True,reversefactor=True,negative_control_gRNAs=self.negative_control_gRNAs)
            else:
                self.size_f=None
        else:
            self.size_f=normalizeCounts(cttab_sel,returnfactor=True,reversefactor=True)
        logging.info('size factor: '+','.join([str(x) for x in self.size_f]))
        #logging.info(self.design_matrix)
        desmat=self.design_matrix
        for (tgid,tginst) in self.allgenedict.iteritems():
            #if tgid=="INO80B":
            tginst.design_mat=desmat
            iteratenbem(tginst,debug=False,estimateeff=False,alpha_val=0.05,size_factor=self.size_f)
            tginst.w_estimate=[]
        deviation(self.allgenedict)

        self.non_PPI_beta_prior_variance=beta_non_PPI_prior_calculation(self.allgenedict)
        pickle.dump(self,open("/%s/%s_self_bayes_init.p" %(self.output_directory,self.output_prefix),'wb'))

    def parameter_reevalution(self):
        self=pickle.load(open("/%s/%s_self_bayes_init.p" %(self.output_directory,self.output_prefix),'rb'))

        logging.info('Reestimate dispersion ...')
        ngenes=0
        for (tgid,tginst) in self.allgenedict.iteritems():
            if ngenes<3000:
                try:
                    #logging.info(tgid)
                    sgrna_wide_dispersion_estimation_MAP_v2(tginst,self.design_matrix)
                    ngenes+=1
                except:
                    pass
        pickle.dump(self,open("/%s/%s_self_para_eval.p" %(self.output_directory,self.output_prefix),'wb'))

    def bayes_fitting(self):
        self=pickle.load(open("/%s/%s_self_para_eval.p" %(self.output_directory,self.output_prefix),'rb'))
        logging.info('Modeling the mean and variance ...')
        if self.genes_varmodeling>0:
            self.mrm=MeanVarModel()
            #self.mrm.model_mean_disp_by_lm(self.allgenedict)
            self.mrm.model_mean_disp_by_glm(self.allgenedict,self.output_prefix,self.size_f)
        else:
            mrm=None
        pickle.dump(self,open("/%s/%s_self_fitting.p" %(self.output_directory,self.output_prefix),'wb'))

    def bayes_major(self):
        self=pickle.load(open("/%s/%s_self_fitting.p" %(self.output_directory,self.output_prefix),'rb'))
        logging.info('Run the algorithm for the second time ...')
        ngene=0
        for (tgid,tginst) in self.allgenedict.iteritems():
            #logging.info(tgid)
            iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=False,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
            ngene+=1

        for (tgid,tginst) in self.allgenedict.iteritems():
            if len(tginst.w_estimate)==0:
                tginst.w_estimate=np.ones(len(tginst.sgrnaid))

        pickle.dump(self,open("/%s/%s_self_bayes_major.p" %(self.output_directory,self.output_prefix),'wb'))

    def constant_optimization(self):
        # PPI
        self=pickle.load(open("/%s/%s_self_bayes_major.p" %(self.output_directory,self.output_prefix),'rb'))
        logging.info('PPI validation...')
        PPI_main(self)
        if self.PPI_diagnosis==True:
            beta_prior_output(self)
        '''
        copy_self=copy.copy(self)
        copy_self.count_table="singlized_{}".format(self.count_table)
        copy_self.output_prefix="{}_mageck_bayes".format(copy_self.count_table)
        copy_self.output_directory="{}/{}".format(os.getcwd(),copy_self.output_prefix)
        if os.path.exists(copy_self.output_directory)==False:
            os.makedirs(copy_self.output_directory)
        if not os.path.isfile("/%s/self_pre_dispersion_fitting_round.p" %copy_self.output_directory):
            os.chdir(self.file_directory)
            sgRNA_production(self.count_table)
            copy_self.allgenedict,copy_self.invalid_gRNA_dict=read_gene_from_file(copy_self.count_table,includesamples=copy_self.include_samples)
            for (tgid,tginst) in copy_self.allgenedict.items():
                tginst.design_mat=copy_self.design_matrix
                iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=False,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
            pickle.dump(copy_self,open("/%s/self_pre_dispersion_fitting_round.p" %copy_self.output_directory,'wb'))
        else:
            copy_self=pickle.load(open("/%s/self_pre_dispersion_fitting_round.p" %copy_self.output_directory,'rb'))
        logging.info("signlized_sgRNA processed!")
        sgRNA_beta=background(copy_self)
        beta1=list(sgRNA_beta.values())
        median=np.median(beta1)
        abs_beta1=[abs(i-median) for i in beta1]
        abs_beta1.sort()
        std_beta1=float(np.percentile(abs_beta1,95))/norm.ppf(0.975)
        for geneid in list(self.allgenedict.keys()):
            sgrnaid=self.allgenedict[geneid].sgrnaid
            temp=[]
            for sgRNA in sgrnaid:
                temp.append(np.log(norm.pdf(sgRNA_beta[sgRNA],median,std_beta1)))
            self.allgenedict[geneid].sgrna_probability=temp
        '''
        logging.info('Estimate selection constant...')
        bayes_selection_constnat_optimization(self)
        pickle.dump(self,open("/%s/%s_self_constant_optimization.p" %(self.output_directory,self.output_prefix),'wb'))

    def bayes_iteration(self):
        logging.info('Iteratioin...')
        PPI_choice=[False,True]

        for PPI in PPI_choice:
            self=pickle.load(open("/%s/%s_self_constant_optimization.p" %(self.output_directory,self.output_prefix),'rb'))
            self.PPI_prior=PPI
            if self.PPI_diagnosis==False and self.PPI_prior==True:
                pass
            else:
                self.outliers_removal=False
                for (tgid,tginst) in self.allgenedict.iteritems():
                    iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=False,PPI_prior=self.PPI_prior,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
                pickle.dump(self,open("/%s/%s_self_bayes_iteration_PPI_%s_outliers_removal_%s.p" %(self.output_directory,self.output_prefix,self.PPI_prior,self.outliers_removal),'wb'))

                #base_index=np.where(~self.design_matrix[:,1:].any(axis=1))[0]
                #transform_list=[1 if i in base_index else -1 for i in range(self.design_matrix.shape[0])]
                #transform_matrix=np.matrix(transform_list)

                for (tgid,tginst) in self.allgenedict.items():
                    bayes_selection(tginst,log_list=self.log_list,selection_constant=self.selection_constant)
                self.outliers_removal=True
                for (tgid,tginst) in self.allgenedict.iteritems():
                    iteratenbem(tginst,debug=False,meanvarmodel=self.mrm,restart=False,PPI_prior=self.PPI_prior,removeoutliers=self.outliers_removal,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
                pickle.dump(self,open("/%s/%s_self_bayes_iteration_PPI_%s_outliers_removal_%s.p" %(self.output_directory,self.output_prefix,self.PPI_prior,self.outliers_removal),'wb'))

    def bayes_output(self):
        os.chdir(self.output_directory)
        for mark in [[True,False],[True,True],[False,False],[False,True],["major","major"]]:
            if mark!=["major","major"]:
                file="%s_self_bayes_iteration_PPI_%s_outliers_removal_%s.p" %(self.output_prefix,mark[0],mark[1])
            else:
                file="%s_self_bayes_major.p" %self.output_prefix
            if file in glob.glob("*"):
                self=pickle.load(open(file,'rb'))
                # permutation
                iteratenbem_permutation(self.allgenedict,nround=self.permutation_round,meanvarmodel=self.mrm,size_factor=self.size_f,beta1_prior_var=self.non_PPI_beta_prior_variance)
                # correct for FDR
                gene_fdr_correction(self.allgenedict,self.adjust_method)

                genefile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.gene_summary.txt'
                sgrnafile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.sgrna_summary.txt'
                logging.info(genefile)
                logging.info(sgrnafile)

                logging.info('Writing gene results to '+genefile)
                logging.info('Writing sgRNA results to '+sgrnafile)
                write_gene_to_file(self.allgenedict,genefile,betalabels=self.beta_labels)
                write_sgrna_to_file(self.allgenedict,sgrnafile)
                plot(genefile)


    def bayes_output_no_permutation(self):
        os.chdir(self.output_directory)
        marks=[['True','False'],['True','True'],['False','False'],['False','True'],["major","major"]]
        #marks=[["major","major"]]
        deviation_output=open("%s_mu_k_deviation.txt" %self.output_prefix,'wb')
        deviation_insert=["PPI","outliers_removal","central_80_percent_mean","median"]
        deviation_output.write("%s\n" %"\t".join(deviation_insert))

        PPI_diagnosis_output=open("%s_mu_k_PPI_diagnosis.txt" %self.output_prefix,'wb')
        PPI_diagnosis_insert=[""]
        PPI_diagnosis_output.write("%s\n" %"\t".join(PPI_diagnosis_insert))
        '''
        os.chdir("%s/doc/" %self.output_directory)
        DE_file=open("%s_DE.txt" %self.count_table)
        DE_file.readline()
        DE_dict=dict()
        for line in DE_file:
            elements=line.strip().split("\t")
            DE_dict[elements[0].upper()]=elements[4]
        '''
        os.chdir(self.output_directory)
        for mark in marks:
            if mark!=["major","major"]:
                file="%s_self_bayes_iteration_PPI_%s_outliers_removal_%s.p" %(self.output_prefix,mark[0],mark[1])
            else:
                file="%s_self_bayes_major.p" %self.output_prefix
            if file in glob.glob("*"):
                self=pickle.load(open(file,'rb'))
                # correct for FDR

                gene_fdr_correction(self.allgenedict,self.adjust_method,wald_only=True)

                genefile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.gene_summary.txt'
                sgrnafile=self.output_prefix+'_PPI_'+str(mark[0])+'_outliers_removal_'+str(mark[1])+'.sgrna_summary.txt'

                logging.info('Writing gene results to '+genefile)
                logging.info('Writing sgRNA results to '+sgrnafile)
                write_gene_to_file(self.allgenedict,genefile,betalabels=self.beta_labels,wald_only=True)
                #write_sgrna_to_file(self.include_samples,self.allgenedict,self.invalid_gRNA_dict,sgrnafile,DE_dict)
                write_sgrna_to_file(self.include_samples,self.allgenedict,self.invalid_gRNA_dict,sgrnafile)
                plot(genefile,negative_control_gRNAs=self.negative_control_gRNAs,wald_only=True)

                deviation_insert=mark+deviation(self.allgenedict)
                deviation_output.write("%s\n" %"\t".join(deviation_insert))

                #PPI_diagnosis_insert=PPI_main(self,final_evaluation=True)
                #PPI_diagnosis_output.write("%s\n" %"\t".join(mark))
                #logging.info(PPI_diagnosis_insert)
                #PPI_diagnosis_output.write("r-square: %s\n" %"\t".join([str(i) for i in PPI_diagnosis_insert[0][0]]))
                #PPI_diagnosis_output.write("regression: %s\n" %"\t".join([str(i) for i in PPI_diagnosis_insert[0][1]]))

def mageck_bayes_main(pvargs=None,parsedargs=None,returndict=False):
    if parsedargs==None:
      initial_args=crisprseq_parseargs()
      args=postargs(initial_args)
    else:
      args=parsedargs
    g=Mageck_Bayes(args)

    g.bayes_init()
    g.parameter_reevalution()

    g.bayes_fitting()
    g.bayes_major()
    g.constant_optimization()

    g.bayes_iteration()

    g.bayes_output_no_permutation()

if __name__ == '__main__':
    try:
        mageck_bayes_main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)
