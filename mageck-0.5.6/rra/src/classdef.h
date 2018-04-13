#ifndef CLASSDEF_H
#define CLASSDEF_H



#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define NDEBUG
#include <assert.h>

#define MAX_NAME_LEN 10000           //maximum length of item name, group name or list name
#define CDF_MAX_ERROR 1E-10        //maximum error in Cumulative Distribution Function estimation in beta statistics
#define MAX_GROUP_NUM 100000       //maximum number of groups
#define MAX_LIST_NUM 1000          //maximum number of list 
#define RAND_PASS_NUM 100          //number of passes in random simulation for computing FDR

#define MAX_WORD_NUM 1000        //maximum number of word

//WL: define boolean variables
//typedef int bool;
//#define true 1
//#define false 0


typedef struct // item definition; i.e., sgRNA
{
	char name[MAX_NAME_LEN];       //name of the item
	int listIndex;                 //index of list storing the item
	double value;                  //value of measurement
	double percentile;             //percentile in the list
	double prob;			//The probability of each sgRNA; added by Wei
  int isChosen;  //whether this sgRNA should be considered in calculation
} ITEM_STRUCT;

typedef struct // group definition; i.e., gene
{
	char name[MAX_NAME_LEN];       //name of the group
	ITEM_STRUCT *items;            //items in the group
	int itemNum;                   //number of items in the group
  int maxItemNum;                // max number of items
	double loValue;                //lo-value in RRA
  double pvalue;                //p value for permutation
	double fdr;                    //false discovery rate
  int isbad;                    //if the lovalue is too low (i.e., higher than the given percentile)
  int goodsgrnas;               //sgRNAs with significant changes
  int controlsgs;               //the numbers of control sgs
} GROUP_STRUCT;

typedef struct //list definition; i.e., gene groups
{
	char name[MAX_NAME_LEN];       //name of the list
	double *values;                //values of items in the list, used for sorting
	int itemNum;                   //number of items in the list
  int maxItemNum;               //max item number
} LIST_STRUCT;


#endif
