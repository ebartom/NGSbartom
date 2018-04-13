//C++ functions
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
using namespace std;

#include "fileio.h"


//split strings
int stringSplit(string str, string delim, vector<string> & v){
  int start=0;
  int pos=str.find_first_of(delim,start);
  v.clear();
  while(pos!=str.npos){
    if(pos!=start)
      v.push_back(str.substr(start,pos-start));
     start=pos+1;
     pos=str.find_first_of(delim,start);
  }
  if(start<str.length())
    v.push_back(str.substr(start,str.length()-start));
  return v.size();
}

int getGroupListNum(char* fileName, GROUP_STRUCT* &groups, LIST_STRUCT* lists, int maxGroupNum, int maxListNum, int &groupNum, int& listNum, map<string,int>& groupNames, map<string,int>& listNames){
  // determine the group number and list number
  // if success, return 0; else, return -1;
  int i,k;
  //char **words, *tmpS, *tmpS2;
  int wordNum;
  int tmpGroupNum, tmpListNum;

  ifstream fh;
  fh.open(fileName);
  if(!fh.is_open()){
    cerr<<"Error opening "<<fileName<<endl;
  }
  string oneline;
  getline(fh,oneline);
  
  vector<string> vwords;
  vector<string> vsubwords;
  wordNum=stringSplit(oneline," \t\r\n\v",vwords);
  //Read the header row to get the sample number
  //fgets(tmpS, MAX_WORD_NUM*(MAX_NAME_LEN+1)*sizeof(char), fh);
  //wordNum = StringToWords(words, tmpS,MAX_WORD_NUM, MAX_NAME_LEN+1,  " \t\r\n\v\f");
  if (wordNum < 4 || wordNum > 6){
    cerr<<"Error: incorrect input file format: <item id> <group id> <list id> <value> [<prob>] [iscounted]\n";
    fh.close();
    return -1;
  }

  //read records of items

  tmpGroupNum = 0;
  tmpListNum = 0;
  //fgets(tmpS, MAX_WORD_NUM*(MAX_NAME_LEN+1)*sizeof(char), fh);
  //wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_NUM, " \t\r\n\v\f");
  //
  getline(fh,oneline);
  wordNum=stringSplit(oneline," \t\r\f\v",vwords);

	
  //WL
  int subWordNum;
  subWordNum=0;

  while ((wordNum>=4 )&&(!fh.eof()))
  {
    
    // separate the group name by ","
    subWordNum=stringSplit(vwords[1],",",vsubwords);
    assert(subWordNum>0);
    
    for(k=0;k<subWordNum;k++){
       //for (i=0;i<tmpGroupNum;i++) if (!strcmp(vsubwords[k].c_str(), groups[i].name))	break;
       string subwstr=vsubwords[k];
       if (groupNames.count(subwstr)==0){
         strcpy(groups[tmpGroupNum].name, vsubwords[k].c_str());
         groups[tmpGroupNum].itemNum = 1;
         groupNames[subwstr]=tmpGroupNum;
         tmpGroupNum ++;
         if (tmpGroupNum >= maxGroupNum){
           printf("Warning: too many groups; maxGroupNum = %d\n. Will try to double the maximum number of groups..", maxGroupNum);
           GROUP_STRUCT* gp2=new GROUP_STRUCT[maxGroupNum*2];
	   memcpy(gp2,groups,maxGroupNum*sizeof(GROUP_STRUCT));
	   delete []groups;
	   maxGroupNum*=2;
	   groups=gp2;

           //return -1;
         }
       }
       else{
          i=groupNames[subwstr];
          groups[i].itemNum++;
       }
    }

    //for (i=0;i<tmpListNum;i++)	if (!strcmp(vwords[2].c_str(), lists[i].name))	break;
    string thisliststr=vwords[2];
    if (listNames.count(thisliststr)==0){
      strcpy(lists[tmpListNum].name, vwords[2].c_str());
      lists[tmpListNum].itemNum = 1;
      listNames[thisliststr]=tmpListNum;
      tmpListNum ++;
      if (tmpListNum >= maxListNum){
        printf("Error: too many lists. maxListNum = %d\n", maxListNum);
        return -1;
      }
    }
    else{
      i=listNames[thisliststr];
      lists[i].itemNum++;
    }


    getline(fh,oneline);
    wordNum=stringSplit(oneline," \t\r\f\v",vwords);
    //fgets(tmpS, MAX_WORD_NUM*(MAX_NAME_LEN+1)*sizeof(char), fh);
    //wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_NUM, " \t\r\n\v\f");
  }
  fh.close();
  
  groupNum=tmpGroupNum;
  listNum=tmpListNum;

  return 0;
}

//Read input file. File Format: <item id> <group id> <list id> <value>. Return 1 if success, -1 if failure
int ReadFile(char *fileName, GROUP_STRUCT* &groups, int maxGroupNum, int *groupNum, 
      LIST_STRUCT *lists, int maxListNum, int *listNum)
{
	//FILE *fh;
	int i,j,k;
	//char **words, *tmpS, *tmpS2;
	int wordNum=0;
	int totalItemNum=0;
	int tmpGroupNum=0, tmpListNum=0;
	//char tmpGroupName[MAX_NAME_LEN], tmpListName[MAX_NAME_LEN], tmpItemName[MAX_NAME_LEN];
	double tmpValue=0;
  string oneline;
  
  //char **subwords;
	double sgrnaProbValue=0.0;
  int sgrnaChosen=1;
  
	ifstream fh;
  vector<string> vwords;
  vector<string> vsubwords;
  int subWordNum=0;
	
  map<string,int> groupNames;
  map<string,int> listNames;
  //assert(tmpS!=NULL && tmpS2!=NULL);
  if(getGroupListNum(fileName, groups, lists, maxGroupNum, maxListNum, tmpGroupNum, tmpListNum, groupNames, listNames)!=0){
    return -1;
  }
  // construct the structure
	for (i=0;i<tmpGroupNum;i++){
    //printf("Group %d items:%d\n",i,groups[i].itemNum);
		groups[i].items = new ITEM_STRUCT[groups[i].itemNum];
		groups[i].maxItemNum = groups[i].itemNum;
		groups[i].itemNum = 0;
	}
	
	for (i=0;i<tmpListNum;i++){
    //printf("List %d items:%d\n",i,lists[i].itemNum);
		lists[i].values = new double[lists[i].itemNum]; 
		lists[i].maxItemNum = lists[i].itemNum;
		lists[i].itemNum = 0;
	}
	
  // actually loading the file
  fh.open(fileName);
  if(!fh.is_open()){
    cerr<<"Error opening "<<fileName<<endl;
  }
	
	//Read the header row to get the sample number
  getline(fh,oneline);
	//fgets(tmpS, MAX_WORD_NUM*(MAX_NAME_LEN+1)*sizeof(char), fh);
	
	//read records of items
	
  getline(fh,oneline);
  
  wordNum=stringSplit(oneline," \t\r\n\v",vwords);
  // fields saved to vwords: sgRNA name, gene name, list name, value
	int skippedsgrna=0;
	while ((wordNum>=4)&&(!fh.eof())){
		//strcpy(tmpListName, words[2]);
		tmpValue = atof(vwords[3].c_str());
		//WL
		//parsing prob column, if available
		if (wordNum > 4){
			sgrnaProbValue=atof(vwords[4].c_str());
    }else{
			sgrnaProbValue=1.0;
		}
		if (wordNum > 5){
			sgrnaChosen=atoi(vwords[5].c_str());
      if(sgrnaChosen==0) skippedsgrna+=1;
    }else{
			sgrnaChosen=1;
		}


		
    //WL: now, search for corresponding list index
    //for (j=0;j<tmpListNum;j++) if (!strcmp(vwords[2].c_str(), lists[j].name))	break;
    string thislistname=vwords[2];
    assert(listNames.count(thislistname)>0);
    j=listNames[thislistname];
    assert(j<tmpListNum);
    // save to list
    if(sgrnaChosen){
		  lists[j].values[lists[j].itemNum] = tmpValue;
 		  lists[j].itemNum ++;
      assert(lists[j].itemNum <=lists[j].maxItemNum);
    }
    
    //WL: search for gene name (group name)
		// WL: separate the group name by ","
    subWordNum=stringSplit(vwords[1],",",vsubwords);
		assert(subWordNum>0);
		//for(k=0;k<subWordNum;k++)
	 for(k=0;k<subWordNum;k++) {
  			//for (i=0;i<tmpGroupNum;i++)	if (!strcmp(vsubwords[k].c_str(), groups[i].name))	break;
  		
  		  assert(groupNames.count(vsubwords[k])>0);
        i=groupNames[vsubwords[k]];
  			assert(i<tmpGroupNum);
  		
        //save group name(gene name)
  			//strcpy(groups[i].items[groups[i].itemNum].name,tmpItemName);
  			strcpy(groups[i].items[groups[i].itemNum].name,vwords[0].c_str());
  			groups[i].items[groups[i].itemNum].value = tmpValue;
  			groups[i].items[groups[i].itemNum].prob= sgrnaProbValue;
  			groups[i].items[groups[i].itemNum].listIndex = j;
  			groups[i].items[groups[i].itemNum].isChosen= sgrnaChosen;
  			groups[i].itemNum ++;
      	assert(groups[i].itemNum <=groups[i].maxItemNum);
  		
    }//end subwordnum 
		
		totalItemNum++;
		//wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, MAX_WORD_NUM, " \t\r\n\v\f");
    getline(fh,oneline);
    wordNum=stringSplit(oneline," \t\r\f\v",vwords);
   	//printf("wordnum:%d,read one line length:%d\n",wordNum,strlen(tmpS));
	}//end loop for fh reading
	
	fh.close();
	//fclose(fh);
	
	printf("Summary: %d sgRNAs, %d genes, %d lists; skipped sgRNAs:%d\n", totalItemNum, tmpGroupNum, tmpListNum,skippedsgrna);
	
	*groupNum = tmpGroupNum;
	*listNum = tmpListNum;
	
	//FreeWords(words, MAX_WORD_NUM);
	//free(tmpS);
	//free(tmpS2);
	//FreeWords(subwords, MAX_WORD_NUM);
	
	return totalItemNum;
	
}

//Save group information to output file. Format <group id> <number of items in the group> <lo-value> <false discovery rate>
int SaveGroupInfo(char *fileName, GROUP_STRUCT *groups, int groupNum)
{
  FILE *fh;
  int i;
  
  fh = (FILE *)fopen(fileName, "w");
  
  if (!fh){
    printf("Cannot open %s.\n", fileName);
    return -1;
  }
  
  fprintf(fh, "group_id\titems_in_group\tlo_value\tp\tFDR\tgoodsgrna\n");
  
  for (i=0;i<groupNum;i++){
    if(groups[i].controlsgs>=groups[i].itemNum){ //skip those that consists of only control sgrnas
      printf("Suppressing the output of gene %s since it is negative ontrol genes.\n",groups[i].name);
      continue;
    }
    if(groups[i].isbad==0){
      fprintf(fh, "%s\t%d\t%10.4e\t%10.4e\t%f\t%d\n", groups[i].name, groups[i].itemNum, groups[i].loValue, groups[i].pvalue,groups[i].fdr,groups[i].goodsgrnas);
    }else{
      fprintf(fh, "%s\t%d\t%10.4e\t%10.4e\t%s\t%d\n", groups[i].name, groups[i].itemNum, groups[i].loValue, groups[i].pvalue,"NA",groups[i].goodsgrnas);
    }
  }

  fclose(fh);

  return 1;
}


