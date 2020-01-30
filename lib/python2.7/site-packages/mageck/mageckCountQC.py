"""MAGeCK Count QC 
   Perform Quality Control test for MAGeCK count module
"""
from __future__ import print_function

import sys
import argparse
import math
import logging
import string

from mageck.testVisualCount import *
from mageck.fileOps import *


def mageckcount_printstat(args,datastat):
  '''
  Write data statistics to PDF file 
  '''
  for (k,v) in datastat.items():
    logging.info('Summary of file '+k+':')
    for (v1,v2) in v.items():
      logging.info(str(v1)+'\t'+str(v2))
  # write to table
  crv=VisualRCount()
  crv.setPrefix(args.output_prefix)
  alllabels={}
  for (fq, fqstat) in datastat.items():
    crv.fastqfile+=[fq]
    if 'label' in fqstat:
      crv.fastqlabels+=[fqstat['label']]
      alllabels[fqstat['label']]=1
    else:
      crv.fastqlabels+=['NA']
    if 'reads' in fqstat:
      crv.reads+=[fqstat['reads']]
    else:
      crv.reads+=[0]
    if 'mappedreads' in fqstat:
      crv.mappedreads+=[fqstat['mappedreads']]
    else:
      crv.mappedreads+=[0]
    if 'totalsgrnas' in fqstat:
      crv.totalsgrnas+=[fqstat['totalsgrnas']]
    else:
      crv.totalsgrnas+=[0]
    if 'zerosgrnas' in fqstat:
      crv.zerocounts+=[fqstat['zerosgrnas']]
    else:
      crv.zerocounts+=[0]
    if 'giniindex' in fqstat:
      crv.gini+=[fqstat['giniindex']]
    else:
      crv.gini+=[0.0]
    if 'NegSelectionQC' in fqstat:
      crv.negselqc+=[fqstat['NegSelectionQC']]
    else:
      crv.negselqc+=[0.0]
    if 'NegSelectionQCPval' in fqstat:
      crv.negselqcpval+=[fqstat['NegSelectionQCPval']]
    else:
      crv.negselqcpval+=[1.0]
    if 'NegSelectionQCPvalPermutation' in fqstat:
      crv.negselqcpvalpermutation+=[fqstat['NegSelectionQCPvalPermutation']]
    else:
      crv.negselqcpvalpermutation+=[1.0]
    if 'NegSelectionQCPvalPermutationFDR' in fqstat:
      crv.negselqcpvalpermutationfdr+=[fqstat['NegSelectionQCPvalPermutationFDR']]
    else:
      crv.negselqcpvalpermutationfdr+=[1.0]
    if 'NegSelectionGene' in fqstat:
      crv.negselqcgene+=[fqstat['NegSelectionGene']]
    else:
      crv.negselqcgene+=[0.0]
  #
  crv.startRTemplate()
  crv.writeCountSummary()
  # write to TXT file
  crv.writeCountSummaryToTxt(args.output_prefix+'.countsummary.txt')
  # write to PDF file
  outcsvfile=args.output_prefix+'.count_normalized.txt'
  crv.insertReadCountBoxPlot(os.path.basename(outcsvfile))
  if len(alllabels) > 1:
    crv.insertPCAPlot(os.path.basename(outcsvfile))
    crv.insertClusteringPlot(os.path.basename(outcsvfile))
  crv.closeRTemplate()
  if hasattr(args,"pdf_report") and args.pdf_report:
    if hasattr(args,"keep_tmp") :
      crv.generatePDF(keeptmp=args.keep_tmp)
    else:
      crv.generatePDF()


def mageckcount_getQC(args,datastat,sgdict):
  """
  Get the degree of negative selection as a QC
  Parameters:
    args
        The parseargs object
    datastat
        A {fastq:{property:value}} object
    sgdict
        A {sequence:(sgID,gene)} object
  """
  # count table
  if args.count_table == None:
    countfile=args.output_prefix+'.count.txt'
  else:
    countfile=args.count_table
  logging.info('Starting the quality control module.')
  # check the pathway used for enrichment
  # GMT file is checked at check_args()
  try:
    with open(args.gmt_file,'r') as f:
      firstline=f.readline().strip().split()
      targetpathway=firstline[0]
      targetpathwayanno=firstline[1]
      pathwaygene=[x.upper() for x in firstline[2:]]
  except:
    logging.error('Cannot open pathway GMT file '+args.gmt_file+' for QC. Please double check your file location and format.')
    return
  # check whether there are overlaps between pathway and library
  genedict={}
  for (k, sv) in sgdict.items():
    genedict[sv[1].upper()]=0
  pathwaygene_sel=[x for x in pathwaygene if x in genedict]
  if len(pathwaygene_sel)>0:
    logging.info('' + str(len(pathwaygene_sel)) + ' out of ' + str(len(pathwaygene)) + ' genes in ' +targetpathway+' are found in the library. These genes will be used for QC.')
  else:
    logging.warning('Found 0 genes in ' + args.gmt_file + ' that apper in your screening library.')
    # return
  
  # determine the control of sample label
  control_label=','.join(args.day0_label)
  #
  logging.info('control label: '+control_label)
  # for each treatment label, run mageck
  for tlabel in args.sample_label:
    if tlabel in args.day0_label:
      continue
    enrichscore=0.0
    enrichscorep=1.0
    enrichscoreppermutation=1.0
    enrichscoreppermutationfdr=1.0
    outputp=args.output_prefix+"_QC."+tlabel+"_vs_"+control_label
    if len(pathwaygene_sel)>0:
      mageckcommand=sys.argv[0]+" test "
      mageckcommand+= " -k " + countfile
      mageckcommand+= " -n " + outputp
      mageckcommand+= " -t " + tlabel + " -c " + control_label
      mageckcommand+= " --norm-method "+args.norm_method
      if args.control_sgrna != None:
        mageckcommand+= " --control-sgrna "+args.control_sgrna
      logging.info('Commands: '+mageckcommand)
      systemcall(mageckcommand,cmsg=False)
      # check for gmt file
      logging.info('Checking for pathway enrichment ...')
      ## using mageck pathway command
      # mageckpathwaycommand=sys.argv[0]+" pathway "
      # mageckpathwaycommand+= " --gene-ranking "+ outputp+".gene_summary.txt "
      # mageckpathwaycommand+= " --gmt-file \""+args.gmt_file+"\""
      # mageckpathwaycommand+= " -n " + outputp
      # mageckpathwaycommand+= " --single-ranking " 
      mageckpathwaycommand="mageckGSEA "
      mageckpathwaycommand+="-c 2 "
      mageckpathwaycommand+="-p 10000 "
      mageckpathwaycommand+= " -g \""+args.gmt_file+"\""
      mageckpathwaycommand+= " -r \""+ outputp+".gene_summary.txt\" "
      mageckpathwaycommand+= " -o \""+ outputp+".pathway_summary.txt\" "
      logging.info('Pathway commands: '+mageckpathwaycommand)
      systemcall(mageckpathwaycommand,cmsg=False)
      # get results from pathway file
      try:
        for line in open(outputp+'.pathway_summary.txt','r'):
          field=line.strip().split()
          if field[0] == targetpathway:
            enrichscore=float(field[2])
            enrichscorep=float(field[3])
            enrichscoreppermutation=float(field[4])
            enrichscoreppermutationfdr=float(field[5])
            break
      except IOError:
        logging.error('Cannot read file generated by MAGeCK pathway module:  '+outputp+'.pathway_summary.txt. For more details, use the --keep-tmp option and track the corresponding log files.')
    # end if
    # save to file
    for (sampleid, sampledict) in datastat.items():
      if sampledict['label'] == tlabel:
        datastat[sampleid]['NegSelectionQC']=enrichscore
        datastat[sampleid]['NegSelectionQCPval']=enrichscorep
        datastat[sampleid]['NegSelectionQCPvalPermutation']=enrichscoreppermutation
        datastat[sampleid]['NegSelectionQCPvalPermutationFDR']=enrichscoreppermutationfdr
        datastat[sampleid]['NegSelectionGene']=len(pathwaygene_sel)
        # logging.info('Setting '+tlabel+' NegSelectionQC to '+str(enrichscore))
    # remove tmp file
    if args.keep_tmp == False:
      systemcall('rm '+outputp+'*;',cmsg=False)



