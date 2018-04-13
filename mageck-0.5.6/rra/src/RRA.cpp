/*
 *  RRA.cpp
 *  Implementation of Robust Rank Aggregation (RRA)
 *
 *  Created by Han Xu and Wei Li.
 *  Modified by Wei Li.
 *  Copyright 2017 Dana Farber Cancer Institute. All rights reserved.
 *
 */
#include "math_api.h"
#include "words.h"
#include "rvgs.h"
#include "rngs.h"
#include "classdef.h"
#include "fileio.h"

//C++ functions
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

int PRINT_DEBUG=0;

const char* RRA_VERSION="0.5.6";

//Global variables

//store control sequences
bool UseControlSeq=false;
map<string,int> ControlSeqMap; //save the control sequence name and their index in ControlSeqPercentile;
double* ControlSeqPercentile;

// used for calculation of lo-values
double* tmpLovarray=NULL;
int nLovarray=-1;

map<string,int> gene_to_skip; //skip gene

double min_percentage_goodsgrna=0.0;
int min_number_goodsgrna=0;

int max_sgnum_permutation_by_group=100; // the maximum number of sgRNA/gene to perform perumtation by group

//Function declarations

//Process groups by computing percentiles for each item and lo-values for each group
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile);

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int start, int end);

//Compute False Discovery Rate based on uniform distribution
int ComputePermutation(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass);
//An improved version by permuting genes from groups
int ComputePermutationbyNumSG(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass);
//Benjamin Hotchberg process to adjust FDR value
int adjustFDR(GROUP_STRUCT *groups, int groupNum);

//print the usage of Command
void PrintCommandUsage(const char *command);


//Compute lo-value based on an array of percentiles
int ComputeLoValue(double *percentiles,     //array of percentiles
   int num,                 //length of array
   double &loValue,         //pointer to the output lo-value
   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
   int &goodsgrna);// # of good sgRNAs

//WL: modification of lo_value computation
int ComputeLoValue_Prob(double *percentiles,     //array of percentiles
   int num,                 //length of array
   double &loValue,         //pointer to the output lo-value
   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
   double *probValue,// probability of each prob, must be equal to the size of percentiles
   int &goodsgrna);

typedef pair<double,int> p_ind;
// Correct FDR
int calculateFDRVector(vector<double> & pin, vector<double>& pout);

//read control sequences
int loadControlSeq(const char* fname){
  ifstream fh;
  fh.open(fname);
  if(!fh.is_open()){
    cerr<<"Error opening "<<fname<<endl;
  }
  string oneline;
  int ncount=0;
  ControlSeqMap.clear();
  while(!fh.eof()){
    getline(fh,oneline);
    ControlSeqMap[oneline]=ncount;
    ncount++;
  }
  fh.close();
  ControlSeqPercentile = new double[ncount];
  for(int i=0;i<ncount;i++) ControlSeqPercentile[i]=-1.0;
  cout<<ControlSeqMap.size()<<" control sequences loaded.\n";
  return 0;
}

/*
 * Main entry
 */
int main (int argc, const char * argv[]) {
  int i,flag;
  GROUP_STRUCT *groups;
  int groupNum;
  LIST_STRUCT *lists;
  int listNum;
  char inputFileName[1000], outputFileName[1000];
  double maxPercentile;
  int rand_passnum;
  bool permutation_by_group=true;
  
  //Parse the command line
  if (argc == 1)
  {
    PrintCommandUsage(argv[0]);
    return 0;
  }
  
  printf("Welcome to RRA v %s.\n", RRA_VERSION);
  inputFileName[0] = 0;
  outputFileName[0] = 0;
  maxPercentile = 0.1;
  rand_passnum=RAND_PASS_NUM;
  min_percentage_goodsgrna=0.0;
  min_number_goodsgrna=0;
  bool has_specified_rand_passnum_parameter=false; // whether users specified --permutation option
  
  for (i=2;i<argc;i++)
  {
     if (strcmp(argv[i-1], "-i")==0){
       strcpy(inputFileName, argv[i]);
     }
     if (strcmp(argv[i-1], "-o")==0){
       strcpy(outputFileName, argv[i]);
     }
     if (strcmp(argv[i-1], "-p")==0){
       maxPercentile = atof(argv[i]);
     }
     if (strcmp(argv[i-1], "--control")==0){
       UseControlSeq=true;
       // load control sequences
       loadControlSeq(argv[i]);
     }
     if (strcmp(argv[i-1], "--permutation")==0){
       rand_passnum= atoi(argv[i]);
       has_specified_rand_passnum_parameter=true;
     }
     if (strcmp(argv[i-1], "--skip-gene")==0){
	string sn(argv[i]);
	cout<<"Skipping gene "<<sn<<" for permutation ..."<<endl;
       gene_to_skip[sn]=0;
     }
     if (strcmp(argv[i-1], "--min-percentage-goodsgrna")==0){
       min_percentage_goodsgrna= atof(argv[i]);
     }
     if (strcmp(argv[i-1], "--min-number-goodsgrna")==0){
       min_number_goodsgrna= atoi(argv[i]);
     }
     if (strcmp(argv[i-1], "--max-sgrnapergene-permutation")==0){
       max_sgnum_permutation_by_group= atoi(argv[i]);
     }
  }
  for (i=1;i<argc;i++)
  {
     if (strcmp(argv[i], "--no-permutation-by-group")==0){
        permutation_by_group=false;
     }
  }
  
  if ((inputFileName[0]==0)||(outputFileName[0]==0))
  {
    cerr<<"Error: input file or output file name not set.\n";
    PrintCommandUsage(argv[0]);
    return -1;
  }
  
  if ((maxPercentile>1.0)||(maxPercentile<0.0))
  {
    cerr<<("Error: maxPercentile should be within 0.0 and 1.0\n");
    return -1;
  }
  
  groups = (GROUP_STRUCT *)malloc(MAX_GROUP_NUM*sizeof(GROUP_STRUCT));
  lists = (LIST_STRUCT *)malloc(MAX_LIST_NUM*sizeof(LIST_STRUCT));
  assert(groups!=NULL);
  assert(lists!=NULL);
  
  printf("Reading input file...\n");
  
  flag = ReadFile(inputFileName, groups, MAX_GROUP_NUM, &groupNum, lists, MAX_LIST_NUM, &listNum);
  
  if (flag<=0){
    cerr<<"\nError: reading ranking file ...\n";
    return -1;
  }
  
  cerr<<("Computing lo-values for each group...\n");
  
  if (ProcessGroups(groups, groupNum, lists, listNum, maxPercentile)<=0)
  {
      cerr<<("\nError: processing groups failed.\n");
      return -1;
  }
  	
  cerr<<("Computing false discovery rate...\n");
  	
  if(rand_passnum*groupNum<100000  && has_specified_rand_passnum_parameter==false){
    rand_passnum=100000/groupNum+1;
    cout<<"Increase the number of permutations to "<<rand_passnum<<" to get precise p values. To avoid this, use the --permutation option.\n";
  }
  
  if(permutation_by_group){
    i=ComputePermutationbyNumSG(groups, groupNum, maxPercentile, rand_passnum*groupNum); 
  }else{
    i=ComputePermutation(groups, groupNum, maxPercentile, rand_passnum*groupNum); 
  }
  if (i<=0)
  {
     cerr<<("\nError: computing FDR failed.\n");
     return -1;
  }
  adjustFDR(groups,groupNum);
  	
  cerr<<("Saving to output file...\n");
  
  if (SaveGroupInfo(outputFileName, groups, groupNum)<=0)
  {
     cerr<<("\nError: saving output file failed.\n");
     return -1;
  }
  	
  cerr<<("RRA completed.\n");
  
  for(i=0;i<groupNum;i++)
    delete[] groups[i].items;
  free(groups);
  
  for (i=0;i<listNum;i++)
  {
     delete[] lists[i].values;
  }
  free(lists);
  if(UseControlSeq){
     delete[] ControlSeqPercentile;
  }
   if(nLovarray>0){
      delete []tmpLovarray;
   }
  
  return 0;
  
}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//print the options of the command
	printf("%s - Robust Rank Aggreation v %s.\n", command, RRA_VERSION);
	printf("usage:\n");
	printf("-i <input data file>. Format: <item id> <group id> <list id> <value> [<probability>] [<chosen>]\n");
	printf("-o <output file>. Format: <group id> <number of items in the group> <lo-value> <false discovery rate>\n");
	printf("-p <maximum percentile>. RRA only consider the items with percentile smaller than this parameter. Default=0.1\n");
	printf("--control <control_sgrna list>. A list of control sgRNA names.\n");
	printf("--permutation <int>. The number of rounds of permutation. Increase this value if the number of genes is small. Default 100.\n");
	printf("--no-permutation-by-group. By default, gene permutation is performed separately, by their number of sgRNAs. Turning this option will perform permutation on all genes together. This makes the program faster, but the p value estimation is accurate only if the number of sgRNAs per gene is approximately the same.\n");
	printf("--skip-gene <gene_name>. Genes to skip from doing permutation. Specify it multiple times if you need to skip more than 1 genes.\n");
	printf("--min-percentage-goodsgrna <min percentage>. Filter genes that have too few percentage of 'good sgrnas', or sgrnas that fall below the -p threshold. Must be a number between 0-1. Default 0 (do not filter genes).\n");
	printf("--min-number-goodsgrna <min number>. Filter genes that have too few number of 'good sgrnas', or sgrnas that fall below the -p threshold. Must be an integer. Default 0 (do not filter genes). \n");
	printf("--max-sgrnapergene-permutation <max number>. Only permute genes by group if the number of sgRNAs per gene is smaller than this number. This will save a lot of time if some regions are targeted by a large number of sgRNAs (usually hundreds). Must be an integer. Default 100. \n");
	printf("example:\n");
	printf("%s -i input.txt -o output.txt -p 0.1 \n", command);
	
}




//Process groups by computing percentiles for each item and lo-values for each group
//groups: genes
//lists: a set of different groups. Comparison will be performed on individual list
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile)
{
  int i,j;
  int listIndex, index1, index2;
  int maxItemPerGroup;
  double *tmpF;
  double *tmpProb;
  bool isallone; // check if all the probs are 1; if yes, do not use accumulation of prob. scores
  
  maxItemPerGroup = 0;

  PRINT_DEBUG=1;
	
  for (i=0;i<groupNum;i++){
    if (groups[i].itemNum>maxItemPerGroup){
       maxItemPerGroup = groups[i].itemNum;
    }
  }
  
  assert(maxItemPerGroup>0);
  
  tmpF = new double[maxItemPerGroup];
  tmpProb = new double [maxItemPerGroup];
  
  for (i=0;i<listNum;i++){
    QuicksortF(lists[i].values, 0, lists[i].itemNum-1);
  }
  
  for (i=0;i<groupNum;i++){
    //Compute percentile for each item
    isallone=true;
    int validsgs=0;
    groups[i].controlsgs=0;
    for (j=0;j<groups[i].itemNum;j++){
      if(groups[i].items[j].isChosen==0) continue;
      listIndex = groups[i].items[j].listIndex;
      index1 = bTreeSearchingF(groups[i].items[j].value-0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);
      index2 = bTreeSearchingF(groups[i].items[j].value+0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);

      groups[i].items[j].percentile = ((double)index1+index2+1)/(lists[listIndex].itemNum*2);
      tmpF[validsgs] = groups[i].items[j].percentile;
      tmpProb[validsgs]=groups[i].items[j].prob;
      if(tmpProb[validsgs]!=1.0){
        isallone=false;
      }
      // save for control sequences
      if(UseControlSeq){
        string sgname(groups[i].items[j].name);
        if(ControlSeqMap.count(sgname)>0){
          int sgindex=ControlSeqMap[sgname];
          ControlSeqPercentile[sgindex]=tmpF[validsgs];
	  groups[i].controlsgs++;
        }
      }//end if
      validsgs++;
    }// end j
    if(validsgs<=1){
      isallone=true;
    }
    //printf("Gene: %s\n",groups[i].name);
    if(isallone){
      ComputeLoValue(tmpF, validsgs, groups[i].loValue, maxPercentile, groups[i].goodsgrnas);
    }else{
      ComputeLoValue_Prob(tmpF, validsgs, groups[i].loValue, maxPercentile, tmpProb,groups[i].goodsgrnas);
    }
    // set skip gene if controlsgs equals the number of sgrnas
    if(UseControlSeq){
      if(groups[i].controlsgs==groups[i].itemNum){
	string sn(groups[i].name);
        gene_to_skip[sn]=1;
	cout<<"Skipping gene "<<sn<<" for permutation ..."<<endl;
      }
    }
  }//end i loop

  delete[] tmpF;
  delete[] tmpProb;
  //check if all control sequences are properly assigned a value
  if(UseControlSeq){
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]<0){
        cerr<<"Warning: sgRNA "<<mit->first<<" not found in the ranked list. \n";
        //ControlSeqPercentile[mit->second]=0.5;
      }
    }
  }//end UseControlSeq
	
	return 1;
}

//Compute lo-value based on an array of percentiles. Return 1 if success, 0 if failure
int ComputeLoValue(double *percentiles,     //array of percentiles
   int num,  //length of array
   double &loValue,   //pointer to the output lo-value
   double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
   int &goodsgrna){   // the number of " good" sgRNAs

  int i;
  double *tmpArray;
  double tmpLoValue, tmpF;
  
  if(num==0){
    loValue=1.00;
    goodsgrna=0;
    return 0;
  }
  if(num>nLovarray){
    delete[] tmpLovarray;
    tmpLovarray=new double[num];
    nLovarray=num;
  }
  //tmpArray = (double *)malloc(num*sizeof(double));
  tmpArray=tmpLovarray;
  if (!tmpArray){
    return -1;
  }

  memcpy(tmpArray, percentiles, num*sizeof(double));

  QuicksortF(tmpArray, 0, num-1);
	
  tmpLoValue = 1.0;
  goodsgrna=0;
	
  for (i=0;i<num;i++){
    // if ((tmpArray[i]>maxPercentile)&&(i>0)) //only calculate the value if at least 1 percentile is smaller than the cutoff
    if ((tmpArray[i]>maxPercentile)){
      if(i==0){}
      else{
        break;
      }
    }else{
      goodsgrna++;
    }
    tmpF = BetaNoncentralCdf((double)(i+1),(double)(num-i),0.0,tmpArray[i],CDF_MAX_ERROR);
    if (tmpF<tmpLoValue){
      tmpLoValue = tmpF;
    }
  }

  loValue = tmpLoValue;

  //free(tmpArray);
	
  return 0;
	
}

//Compute lo-value based on an array of percentiles, by considering the probability of each sgRNAs. 
//Return 1 if success, 0 if failure
//Modified by Wei Li
int ComputeLoValue_Prob(double *percentiles,     //array of percentiles
  int num,                 //length of array
  double &loValue,         //pointer to the output lo-value
  double maxPercentile,   //maximum percentile, computation stops when maximum percentile is reached
  double *probValue,
  int &goodsgrna) {// probability of each prob, must be equal to the size of percentiles

  int i, pid;
  double *tmpArray;
  double tmpLoValue, tmpF;

  // used to calculate cumulative probability
  int real_i;
  int c_num;
  double c_prob;
  double accuLoValue;

  if(num==0){
    loValue=1.00;
    goodsgrna=0;
    return 0;
  }

  if(num>nLovarray){
    delete[] tmpLovarray;
    tmpLovarray=new double[num];
    nLovarray=num;
  }
  //tmpArray = (double *)malloc(num*sizeof(double));
  tmpArray=tmpLovarray;

  if (!tmpArray){
    return -1;
  }

  memcpy(tmpArray, percentiles, num*sizeof(double));

  QuicksortF(tmpArray, 0, num-1);
  
  goodsgrna=0;
  for (i=0;i<num;i++){
    if ((tmpArray[i]>maxPercentile)){
      if(i==0){}else{
        break;
      }
    }else{
      goodsgrna++;
    }
  }

  
  if(PRINT_DEBUG) printf("probs:");
  for(i=0;i<num;i++)
  {
    if(PRINT_DEBUG) printf("%f,",probValue[i]);
  }
  if(PRINT_DEBUG) printf("\n");

  accuLoValue=0.0;
  for (pid=1;pid<(1<<num);pid++)
  {
    // decoding the selection
    tmpLoValue = 1.0;
    c_num=0;
    real_i=0;
    c_prob=1.0;
    for(i=0;i<num;i++)
    {
      if ( (1<<i)&pid )
      {
        c_num = c_num +1;
        c_prob= c_prob*probValue[i];
      }
      else
      {
        c_prob=c_prob*(1.0-probValue[i]);
      }
    }
    for (i=0;i<num;i++)
    {
      if ( ((1<<i)&pid) ==0)
      {// if this sgRNA is selected?
        continue;
      }
      if ((tmpArray[i]>maxPercentile)&&(i>0))
      {
        break;
      }
      // Beta (a,b, 0, frac)
      tmpF = BetaNoncentralCdf((double)(real_i+1),(double)(c_num-i),0.0,tmpArray[i],CDF_MAX_ERROR);
      if (tmpF<tmpLoValue)
      {
        tmpLoValue = tmpF;
      }
      real_i = real_i + 1;
     }
    if(PRINT_DEBUG) printf("pid: %d, prob:%e, score: %f\n",pid,c_prob,tmpLoValue);
    accuLoValue=accuLoValue+tmpLoValue*c_prob;
  }
  if(PRINT_DEBUG) printf("total: %f\n",accuLoValue);

  loValue = accuLoValue;

  //free(tmpArray);
  return 0;
	
}

//Benjamin Hotchberg process to adjust FDR value
int adjustFDR(GROUP_STRUCT *groups, int groupNum){
  int i=0;
  //FDR calcuation
  int ngn=0;
  for(int i=0;i<groupNum;i++){
    if(groups[i].goodsgrnas*1.0/groups[i].itemNum < min_percentage_goodsgrna 
		    || groups[i].goodsgrnas<min_number_goodsgrna ){
      groups[i].isbad=1;
    }else{
      ngn++;
    }
  }
  cout<<"Number of genes under FDR adjustment: "<<ngn<<endl;
  vector<double> pvalue_in;
  for(i=0;i<groupNum;i++){
    if(groups[i].isbad==0){
      pvalue_in.push_back(groups[i].pvalue);
    }
  }
  vector<double> pvalue_out(pvalue_in.size());
  
  calculateFDRVector(pvalue_in,pvalue_out);
  int ni=0;
  for(i=0;i<groupNum;i++){
    if(groups[i].isbad==0){
      groups[i].fdr=pvalue_out[ni];
      ni++;
    }else{
      groups[i].fdr=1.0;
    }
  }
  return 0;
}

//Compute False Discovery Rate based on uniform distribution
//All genes with different #gRNAs are permuted differently
int ComputePermutationbyNumSG(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass)
{
  int i,j,k;
  double *tmpPercentile;
  int maxItemNum = 0;
  int scanPass = numOfRandPass/groupNum+1;
  double *randLoValue;
  int randLoValueNum;
  
  //WL
  double *tmpProb;

  map<int,int> itemNumMap; // a map structure to collect different # sgRNAs
	
  for (i=0;i<groupNum;i++){
    groups[i].isbad=0;
    string gn(groups[i].name);
    groups[i].pvalue=1.0;
    groups[i].fdr=1.0;

    if (gene_to_skip.count(gn)>0) continue;
    if(itemNumMap.count(groups[i].itemNum)==0){
      itemNumMap[groups[i].itemNum]=0;
    }
    itemNumMap[groups[i].itemNum]++;
    if (groups[i].itemNum>maxItemNum){
      maxItemNum = groups[i].itemNum;
    }
  }

  assert(maxItemNum>0);
	
  //tmpPercentile = (double *)malloc(maxItemNum*sizeof(double));
  //tmpProb= (double *)malloc(maxItemNum*sizeof(double));
  tmpPercentile=new double[maxItemNum];
  tmpProb=new double[maxItemNum];	

  randLoValueNum = groupNum*scanPass;

  assert(randLoValueNum>0);

  //randLoValue = (double *)malloc(randLoValueNum*sizeof(double));
  randLoValue=new double[randLoValueNum];
	

  PlantSeeds(123456);
  
  PRINT_DEBUG=0;
  
  // set up control sequences
  int n_control=0;
  double* control_prob_array=NULL;
  double ufvalue=0.5;
  int rand_ctl_index=0;
  int tmp_int;
  if(UseControlSeq){
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        n_control++;
      }
    }
    control_prob_array=new double[n_control];
    n_control=0;
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        control_prob_array[n_control]=ControlSeqPercentile[mit->second];
        n_control++;
      }
    }
    cout<<"Total # control sgRNAs: "<<n_control<<endl;
  }
  
  for(map<int,int>::iterator mii=itemNumMap.begin();mii!=itemNumMap.end();mii++){
    int sgrnanum=mii->first;
    if(sgrnanum<=max_sgnum_permutation_by_group){// to avoid permutations on genes with too many sgRNAs (take too much time)
      randLoValueNum = 0;
      cout<<"Permuting genes with "<<sgrnanum<<" sgRNAs..."<<endl;
      for (i=0;i<scanPass;i++){
        for (j=0;j<groupNum;j++){
          int validsgs=0;
          for (k=0;k<sgrnanum;k++)
          {
            ufvalue=Uniform(0.0, 1.0);
            if(UseControlSeq){
              rand_ctl_index=(int)(n_control*ufvalue);
              if(rand_ctl_index>=n_control) rand_ctl_index=n_control-1;
              tmpPercentile[validsgs]=control_prob_array[rand_ctl_index];
            }else{
              tmpPercentile[validsgs] = ufvalue;
            }
            //tmpProb[validsgs]=1.0; // groups[j].items[k].prob;
            //if(tmpProb[validsgs]!=1.0){}
            validsgs++;
          } //end for k
            		
          ComputeLoValue(tmpPercentile, validsgs,randLoValue[randLoValueNum], maxPercentile, tmp_int);
          randLoValueNum++;
        }// end for j
      }//end for i
      QuicksortF(randLoValue, 0, randLoValueNum-1);
    }// end if
    // assigning p value
    //
    for (i=0;i<groupNum;i++){
      if(groups[i].itemNum==sgrnanum){
        string gn(groups[i].name);
        if (gene_to_skip.count(gn)>0) continue;
        if(groups[i].isbad==0){
          groups[i].pvalue=(double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
             +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
             /2/randLoValueNum;
        }else{
          groups[i].pvalue=1.0;
        }
      }
    }

  }//end mii iteration
	
  QuickSortGroupByLoValue(groups, 0, groupNum-1);

  
  delete []tmpPercentile;
  delete []tmpProb;
  delete []randLoValue;
  
  if(UseControlSeq){
    delete[] control_prob_array;
  }

  return 1;
}

//correct the FDR
int calculateFDRVector(vector<double> & pin, vector<double>& pout){
  pout.clear();
  pout.resize(pin.size());
  if(pin.size()==0) return 0;
  vector<p_ind> sortp(pin.size());
  for(int i=0;i<pin.size();i++){
    sortp[i]=p_ind(pin[i],i);
  }
  sort(sortp.begin(),sortp.end());
  for(int i=0;i<sortp.size();i++){
    sortp[i].first=sortp[i].first/(i+1.0)*((int)pin.size());
  }
  if(sortp[sortp.size()-1].first>1.0){
    sortp[sortp.size()-1].first=1.0;
  }
  for(int i=sortp.size()-2;i>=0;i--){
    if(sortp[i].first>sortp[i+1].first){
      sortp[i].first=sortp[i+1].first;
    }
  }
  for(int i=0;i<sortp.size();i++){
    pout[sortp[i].second]=sortp[i].first;
  }
  return 0;
}



//Compute False Discovery Rate based on uniform distribution
//All genes with different #gRNAs are permuted
int ComputePermutation(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass)
{
  int i,j,k;
  double *tmpPercentile;
  int maxItemNum = 0;
  int scanPass = numOfRandPass/groupNum+1;
  double *randLoValue;
  int randLoValueNum;
  
  //WL
  double *tmpProb;
  double isallone;
	
  for (i=0;i<groupNum;i++){
    if (groups[i].itemNum>maxItemNum){
      maxItemNum = groups[i].itemNum;
    }
    groups[i].isbad=0;
  }

  assert(maxItemNum>0);
	
  //tmpPercentile = (double *)malloc(maxItemNum*sizeof(double));
  //tmpProb= (double *)malloc(maxItemNum*sizeof(double));
  tmpPercentile=new double[maxItemNum];
  tmpProb=new double[maxItemNum];	

  randLoValueNum = groupNum*scanPass;

  assert(randLoValueNum>0);

  //randLoValue = (double *)malloc(randLoValueNum*sizeof(double));
  randLoValue=new double[randLoValueNum];
	
  randLoValueNum = 0;

  PlantSeeds(123456);
  
  PRINT_DEBUG=0;
  
  // set up control sequences
  int n_control=0;
  double* control_prob_array=NULL;
  double ufvalue=0.5;
  int rand_ctl_index=0;
  int tmp_int;
  if(UseControlSeq){
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        n_control++;
      }
    }
    control_prob_array=new double[n_control];
    n_control=0;
    for(map<string,int>::iterator mit = ControlSeqMap.begin(); mit != ControlSeqMap.end(); mit++){
      if(ControlSeqPercentile[mit->second]>=0){
        control_prob_array[n_control]=ControlSeqPercentile[mit->second];
        n_control++;
      }
    }
    cout<<"Total # control sgRNAs: "<<n_control<<endl;
  }
  
  for (i=0;i<scanPass;i++){
    for (j=0;j<groupNum;j++){
      isallone=true;
      int validsgs=0;
      for (k=0;k<groups[j].itemNum;k++)
      {
        if(groups[j].items[k].isChosen==0) continue;
        ufvalue=Uniform(0.0, 1.0);
        if(UseControlSeq){
          rand_ctl_index=(int)(n_control*ufvalue);
          if(rand_ctl_index>=n_control) rand_ctl_index=n_control-1;
          tmpPercentile[validsgs]=control_prob_array[rand_ctl_index];
        }else{
          tmpPercentile[validsgs] = ufvalue;
        }
        tmpProb[validsgs]=groups[j].items[k].prob;
        if(tmpProb[validsgs]!=1.0)
        {
          isallone=false;
        }
        validsgs++;
      } //end for k
      if(validsgs<=1)
        isallone=true;
			
      if(isallone){
        ComputeLoValue(tmpPercentile, validsgs,randLoValue[randLoValueNum], maxPercentile, tmp_int);
      }else{
        ComputeLoValue_Prob(tmpPercentile, validsgs,randLoValue[randLoValueNum], maxPercentile,tmpProb,tmp_int);
      }
      randLoValueNum++;
    }// end for j
  }//end for i
	
  QuicksortF(randLoValue, 0, randLoValueNum-1);
  QuickSortGroupByLoValue(groups, 0, groupNum-1);
  
  for (i=0;i<groupNum;i++){
    //groups[i].fdr = (double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
    // +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
    // /2/randLoValueNum/((double)i+0.5)*groupNum;
    if(groups[i].isbad==0){
      groups[i].pvalue=(double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
         +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
         /2/randLoValueNum;
      // groups[i].fdr = groups[i].pvalue/((double)i+1.0)*goodGroupNum;
      // indexval[goodindex]=i;
      // goodindex++;
    }else{
      groups[i].pvalue=1.0;
      // groups[i].fdr=1.0;
    }
  }
	
 

  //free(tmpPercentile);
  //free(tmpProb);
  //free(randLoValue);
  delete []tmpPercentile;
  delete []tmpProb;
  delete []randLoValue;
  
  if(UseControlSeq){
    delete[] control_prob_array;
  }

  return 1;
}

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int lo, int hi)
{
  int i=lo, j=hi;
  GROUP_STRUCT tmpGroup;
  double x=groups[(lo+hi)/2].loValue;

  if (hi<lo)
  {
    return;
  }

  //  partition
  while (i<=j)
  {    
    while ((groups[i].loValue<x)&&(i<=j))
    {
      i++;
    }
    while ((groups[j].loValue>x)&&(i<=j))
    {
      j--;
    }
    if (i<=j)
    {
       memcpy(&tmpGroup,groups+i,sizeof(GROUP_STRUCT));
       memcpy(groups+i,groups+j,sizeof(GROUP_STRUCT));
       memcpy(groups+j,&tmpGroup,sizeof(GROUP_STRUCT));
       i++; j--;
     }
   } 

    //  recursion
    if (lo<j) QuickSortGroupByLoValue(groups, lo, j);
    if (i<hi) QuickSortGroupByLoValue(groups, i, hi);

}
