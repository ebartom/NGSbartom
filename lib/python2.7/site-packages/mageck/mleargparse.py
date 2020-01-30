'''
Argument parsing of the MLE approach

Warning: the mageckmle_parseargs() function is only used in a standalone MLE program. For mageck mle sub-command, the argument parsing is done through argsParser.py.
'''
from __future__ import print_function
import sys
import argparse
import logging


def mageckmle_postargs(args):
  '''
  post-processing of argument parsing
  '''
   # configure logging information
  
  logging.basicConfig(level=10,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    # stream=sys.stderr,
    filename=args.output_prefix+'.log',
    filemode='w'
  )
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  # set a format which is simpler for console use
  formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
  #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
  # tell the handler to use this format
  console.setFormatter(formatter)
  # add the handler to the root logger
  logging.getLogger('').addHandler(console)
  
  logging.info('Parameters: '+' '.join(sys.argv))
 
  from mageck.mledesignmat import parse_designmat,parse_designmat_from_day0
  
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
  if args.design_matrix != None:
    (desmat,sampleid,betalabel)=parse_designmat(args.design_matrix)
  else:
    (desmat,sampleid,betalabel)=parse_designmat_from_day0(args)
  #logging.info(','.join(sampleid))
  #logging.info(','.join(betalabel))
  #print(desmat)
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
  #
  # log design matrix and column, row labels
  logging.info('Design matrix:')
  for desmat_1line in str(desmat).split('\n'):
    logging.info(desmat_1line)
  if args.beta_labels != None:
    logging.info('Beta labels:'+','.join(args.beta_labels))
  if args.include_samples != None:
    logging.info('Included samples:'+','.join(args.include_samples))
  
  return args


def mageckmle_parseargs(pvargs=None):
  '''
  Main function of argument parsing for MAGeCK MLE
  '''

  parser=argparse.ArgumentParser(description='mageck MLE: a maximum likelihood approach for analyzing genome-wide CRISPR/Cas9 knockout screening data.')
  
  # optional parameters
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')
  parser.add_argument('--genes-varmodeling',type=int,default=2000,help='The number of genes for mean-variance modeling. Default 1000.')
  parser.add_argument('--permutation-round',type=int,default=10,help='The rounds for permutation (interger). The permutation time is (# genes)*x for x rounds of permutation. Suggested value: 100 (may take longer time). Default 10.')
  parser.add_argument('-i', '--include-samples', help='Specify the sample labels if the design matrix is not given by file in the --design-matrix option. Sample labels are separated by ",", and must match the labels in the count table.')
  parser.add_argument('-b', '--beta-labels', help='Specify the labels of the variables (i.e., beta), if the design matrix is not given by file in the --design-matrix option. Should be separated by ",", and the number of labels must equal to (# columns of design matrix), including baseline labels. Default value: "bata_0,beta_1,beta_2,...".')
  parser.add_argument('--remove-outliers', action='store_true', help='Try to remove outliers. Turning this option on will slow the algorithm.')
  parser.add_argument('--adjust-method',choices=['fdr','holm','pounds'],default='fdr',help='Method for sgrna-level p-value adjustment, including false discovery rate (fdr) or holm\'s method (holm). Default fdr. (for gene level, always use fdr.)')
  parser.add_argument('--sgrna-efficiency',help='An optional file of sgRNA efficiency prediction. The efficiency prediction will be used as an initial guess of the probability an sgRNA is efficient. Must contain at least two columns, one containing sgRNA ID, the other containing sgRNA efficiency prediction.')
  parser.add_argument('--sgrna-eff-name-column',type=int,default=0,help='The sgRNA ID column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 0 (the first column).')
  parser.add_argument('--sgrna-eff-score-column',type=int,default=1,help='The sgRNA efficiency prediction column in sgRNA efficiency prediction file (specified by the --sgrna-efficiency option). Default is 1 (the second column).')
  parser.add_argument('--update-efficiency',action='store_true',help='Iteratively update sgRNA efficiency during EM iteration.')
  # required parameters
  parser.add_argument('-k','--count-table',required=True,help='Provide a tab-separated count table. Each line in the table should include sgRNA name (1st column), target gene (2nd column) and read counts in each sample.')
  parser.add_argument('-d','--design-matrix',required=True,help='Provide a design matrix, either a file name or a quoted string of the design matrix. For example, "1,1;1,0". The row of the design matrix must match the order of the samples in the count table (if --include-samples is not specified), or the order of the samples by the --include-samples option.')
  #
  if pvargs == None:
    args=parser.parse_args()
  else:
    args=parser.parse_args(pvargs)
  return args
  

