#include <iostream>
#include <fstream>
#include <string>
#include "tclap/CmdLine.h"
#include <map>
#include <vector>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <cstdlib>
#include "math_api.h"

using namespace std;
using namespace TCLAP;

struct CallingArgs{
  string gmtfile;
  string rankfile;
  string outputfile;
  string pname;
  int perms;
  bool needsort;
  int scorecolumn; //column for scoring
  bool reverse; // whether we need to reverse the order of the gene 
};

typedef map<string, vector<string> > pathmap;

//typedef pair<float,int> gs_ret;
struct gs_ret{
  float gscore; // Enrichment score
  int gindex; // index; how many genes are below the cutoff?
  double p_gsea; // p value calculated by GSEA
  double pval; // p value calculated by permutation
  double fdr; // FDR
  string name;
  
  gs_ret(float g, int i): gscore(g), gindex(i){
    p_gsea=1.0;
    pval=1.0;
    fdr=1.0;
    name="";
  }
};

//typedef pair<gs_ret,string> pds;

//typedef pair<double,pds>  d_pds; // a (pval,((score,index),name)) structure

struct dpds_less_key{
  //inline bool operator() (const d_pds& d1, const d_pds& d2){
  //  return (d1.first < d2.first) || ( d1.first == d2.first && d1.second > d2.second );
  //}
  inline bool operator() (const gs_ret& d1, const gs_ret& d2){
    return (d1.pval < d2.pval) || ( d1.pval== d2.pval && d1.gscore> d2.gscore);
  }
};
//map storing gene name and their index
map<string,int> GNAMEINDEX;
//pathways
pathmap PM;
vector<string> GNAME;
vector<double> GSCORE;
//vector<int> INDV;
int* INDV;


const string GSEA_VERSION="0.5.6";

int parseArguments(int argc, char* argv[],CallingArgs& args){
  try{
    CmdLine cmd("mageckGSEA: A fast implementation of GSEA enrichment test.",' ',GSEA_VERSION);
    ValueArg<string> gmtfilearg("g","gmt_file","The pathway annotation in GMT format.",true,"","gmt_file");
    cmd.add(gmtfilearg);
    ValueArg<string> rankfilearg("r","rank_file","Rank file. The first column of the rank file must be the gene name.",true,"","rank_file");
    cmd.add(rankfilearg);
    ValueArg<string> outputfilearg("o","output_file","The name of the output file. Use - to print to standard output.",false,"-","output_file");
    cmd.add(outputfilearg);
    ValueArg<string> pnamearg("n","pathway_name","Name of the pathway to be tested. If not found, will test all pathways.",false,"","pathway_name");
    cmd.add(pnamearg);
    ValueArg<int> permarg("p","perm_time","Permutations, default 1000.",false,1000,"perm_time");
    cmd.add(permarg);
    ValueArg<int> scorecolumnarg("c","score_column","The column for gene scores. If you just want to use the ranking of the gene (located at the 1st column), use 0. The column number starts from 0. Default: 0.",false,0,"score_column");
    cmd.add(scorecolumnarg);
    SwitchArg hassarg("s","sort_byp","Sort the pathways by p value.");
    cmd.add(hassarg);
    SwitchArg reversearg("e","reverse_value","Reverse the order of the gene.");
    cmd.add(reversearg);
    
    cmd.parse(argc,argv);
    
    args.gmtfile=gmtfilearg.getValue();
    args.rankfile=rankfilearg.getValue();
    args.outputfile=outputfilearg.getValue();
    args.pname=pnamearg.getValue();
    args.perms=permarg.getValue();
    args.needsort=hassarg.getValue();
    args.scorecolumn=scorecolumnarg.getValue();
    args.reverse=reversearg.getValue();
    
  }catch(ArgException &e){
    cerr<<"error: "<<e.error()<<" for arg "<<e.argId();
    return -1;
  }
  return 0;
}

int parsegmt(string gmtfile,pathmap& pm){
  ifstream ifs(gmtfile.c_str());
  if(!ifs.is_open()){  
    cerr<<"Error opening mutation file "<<gmtfile<<endl;
    return -1;
  }
  string oneline;
  int ncount=0;
  while(true){
    getline(ifs,oneline);
    ncount++;
    if( ifs.eof())break;
    stringstream ss(oneline);
    int nfield=0;
    string ssf;
    string pname;
    vector<string> vs;
    while(ss.good()){
      nfield++;
      ss>>ssf;
      if(nfield==1){
        pname=ssf;
      }
      if(nfield>2){
        vs.push_back(ssf);
      }
    }
    pm[pname]=vs;
  }
  ifs.close();
  return 0;
}

int parserankfile(string rankf,vector<string>& gname, vector<double>& gscore, CallingArgs &args){

  double providecor=false; //whether the input file provides scores for ranking
  if(args.scorecolumn>0){
    providecor=true;
  }
  int ncount=0;
  
  ifstream ifs(rankf.c_str());
  if(!ifs.is_open()){  
    cerr<<"Error opening rank file "<<rankf<<endl;
    return -1;
  }
  string oneline;
  while(true){
    getline(ifs,oneline);
    ncount++;
    if( ifs.eof())break;
    stringstream ss(oneline);
    int nfield=0;
    string ssf;
    string pname;
    double corval=0.0;
    while(ss.good()){
      nfield++;
      ss>>ssf;
      if(nfield==1){
        pname=ssf;
      }else{
        if(providecor){
          if(nfield==args.scorecolumn+1){
            corval=atof(ssf.c_str());
          }
        }else{
          corval=ncount;
        }
      }
    }
    if(nfield<args.scorecolumn+1){
	cerr<<"Error: cannot find values corresponding to column "<<args.scorecolumn<<" in line "<<ncount<<". Double check your input format."<<endl;
        exit(-1);
    }
    //rmap[pname]=corval;
    gname.push_back(pname);
    gscore.push_back(corval);
  }
  ifs.close();
  //Normalize the values
  if(!providecor){
    for(int i=0;i<gscore.size();i++){
      double t=gscore[i]*2.0/ncount-1.0;
      gscore[i]=t>0?t:-1*t;//abs
    }
  }else{
    typedef pair<double,string> mp;
    vector<mp> gscoreid(gscore.size());
    for(int i=0;i<gscore.size();i++){
	    gscoreid[i]=mp(gscore[i],gname[i]);
    }
    sort(gscoreid.begin(),gscoreid.end());
    for(int i=0;i<gscoreid.size();i++){
      double t=i*2.0/ncount-1.0;
      gscore[i]=t>0?t:-1*t; //abs
      gname[i]=gscoreid[i].second;
      //gscore[gscoreid[i].second]=-1*t; 
    }
  }
  if(args.reverse){
	reverse(gname.begin(),gname.end());
  }
  return 0;
}

//double getgseascore(vector<double> & gscore, vector<int>& ind){
double getgseascore(vector<double> & gscore, int* ind){
  double nf=0.0;
  int n=gscore.size();
  double sumcrv=0.0;
  for(int i=0;i<n;i++){
    if(ind[i]==1){
      nf+=1.0;
      sumcrv+=gscore[i];
    }
  }
  cout<<"Sum:"<<sumcrv<<endl;
  double decvalue=-1.0/(n-nf);
  double gsvalue=0;
  double gsmax=-1000000.0;
  for(int i=0;i<n;i++){
    double incval=gscore[i]/sumcrv;
    if(ind[i]==1){
      gsvalue+=incval;
    }else{
      gsvalue+=decvalue;
    }
    if(gsvalue>gsmax) gsmax=gsvalue;
  }
  
  return gsmax;
}

/*
 * Get the GSEA score
 * Parameters:
 * 	gscore
 * 		The score of each gene
 * 	ind
 * 		The index of the gene in the pathway
 * Return value:
 * 	A gs_ret structure that includes (gsmax,gsmax_ind) value, where
 * 	gsmax
 * 		The maximum of the calculated GSEA score, which is the GSEA score
 * 	gsmax_ind
 * 		(Return value) the index of the maximum GSEA score 
 */
gs_ret getgseascore_fast(vector<double> & gscore, vector<int> & ind){
  double nf=(double)ind.size();
  int n=gscore.size();
  double sumcrv=0.0;
  for(int i=0;i<ind.size();i++){
      //cout<<i<<":"<<ind[i]<<":"<<GNAME[ind[i]]<<",";
      sumcrv+=gscore[ind[i]];
  }
  //cout<<"\nSum:"<<sumcrv<<endl;
  double decvalue=-1.0/(n-nf);
  double gsvalue=0;
  double gsmax=-1000000.0;
  int gsmax_ind=0; // the index where the maximum value achieves
  double sumincval=0.0;
  for(int i=0;i<ind.size();i++){
    int ni=ind[i];
    double incval=gscore[ni]/sumcrv;
    sumincval+=incval;
    //if(ind[i]==1){
    //  gsvalue+=incval;
    //}else{
    //  gsvalue+=decvalue;
    //}
    gsvalue=decvalue*(ni- (i)) + sumincval; 
    if(gsvalue>gsmax) {
      gsmax=gsvalue;
      gsmax_ind=i;
    }
  }
  
  return gs_ret(gsmax,gsmax_ind);
}



gs_ret getscoreandp( string pathname, double & pval, CallingArgs & args){

  int nsz=GNAME.size();
  if(PM.count(pathname)<0){
    cerr<<"Warning: "<<pathname<<" not found.\n";
  }
  double tgscore=0;
  // old method
  //for(int i=0;i<nsz;i++) INDV[i]=0;
  //for(int i =0; i< PM[pathname].size();i++){
  //  string cgname=PM[pathname][i];
  //  INDV[GNAMEINDEX[cgname]]=1;
  //}
  //tgscore=getgseascore(GSCORE,INDV);
  vector<int> index;
  map<int,int> checkred;
  for(int i =0; i< PM[pathname].size();i++){
    string cgname=PM[pathname][i];
    int gidx=GNAMEINDEX[cgname];
    if(checkred.count(gidx)==0){
      //cout<<cgname<<":"<<GNAMEINDEX[cgname]<<endl;
      index.push_back(GNAMEINDEX[cgname]);
      checkred[gidx]=0;
    }
  }
  sort(index.begin(),index.end());
  gs_ret gr=getgseascore_fast(GSCORE,index);
  tgscore=gr.gscore;
  
  
  // get the RRA score
  int num=index.size();
  double listn=GSCORE.size();
  double percentile=index[gr.gindex]*1.0/listn;
  //cout<<gr.gindex<<","<<num-gr.gindex<<":"<<index[gr.gindex]<<"/"<<listn<<"="<<percentile<<endl;
  gr.p_gsea= BetaNoncentralCdf((double)(gr.gindex+1),(double)(num-gr.gindex),0.0,percentile,1.0e-10);

   // Permutation
  int ngt=0;
  for(int p=0;p<args.perms;p++){
    // OLD mehod
    //for(int i=0;i<nsz;i++)INDV[i]=0;
    //int nsp=0;
    //while(nsp<PM[pathname].size()){
    //  int cid=int(rand()*1.0/RAND_MAX*nsz);
    //  if(INDV[cid]==0){
    //    INDV[cid]=1;
    //    nsp++;
    //  }
    //}
    //double ptgscore=getgseascore(GSCORE,INDV);
    
    // NEW method
    checkred.clear();
    int nsp=0;
    while(nsp<index.size()){
      int cid=int(rand()*1.0/RAND_MAX*nsz);
      if(checkred.count(cid)==0){
        checkred[cid]=1;
        index[nsp]=cid;
        nsp++;
      }
    }
    sort(index.begin(),index.end());
    gs_ret gr2=getgseascore_fast(GSCORE,index);
    double ptgscore=gr2.gscore;
    if(ptgscore>tgscore) ngt++;
  }
  //cout<<"Permutation: "<<ngt<<"/"<<args.perms<<" ("<<ngt*1.0/args.perms<<")"<<endl;
  pval=ngt*1.0/args.perms;
  gr.pval=pval;

  return gr;
}



int init_vec(CallingArgs &args){
  cerr<<"GMT file: "<<args.gmtfile<<endl;
  parsegmt(args.gmtfile,PM);
  cerr<<"Rank file: "<<args.rankfile<<endl;
  parserankfile(args.rankfile,GNAME,GSCORE,args);
  cerr<<"Output file: "<<args.outputfile<<endl;
  cerr<<"pathway size:"<<PM.size()<<endl;
  
  cerr<<"rank size:"<<GNAME.size()<<endl;
  cerr<<"Permutations:"<<args.perms<<endl;
  
  srand(time(NULL));
  // create the index
  for(int i=0;i<GNAME.size();i++)
    GNAMEINDEX[GNAME[i]]=i;

  //INDV=vector<int>(GNAME.size(),0);
  INDV=new int[GNAME.size()];
  return 0;
}

int end_vec(CallingArgs &args){
  delete []INDV;
  return 0;
}



// test all pathways in GSEA
int testallgseas(CallingArgs & args){
  //
  double tgscore=0;
  double tgp=0.0;
  vector<gs_ret> vp; // d_pds is a (pval,((score,index),name)) structure

  ofstream ofs;
  if(args.outputfile!="-"){
    ofs.open(args.outputfile.c_str());
    if(!ofs.is_open()){  
      cerr<<"Error opening output file "<<args.outputfile<<endl;
      return -1;
    }
  }
  ostream &gos=(args.outputfile=="-")?cout:ofs;

  gos<<"Pathway\tSize\tES\tp\tp_permutation\tFDR\tRanking\tHits\tLFC\n";
  for(pathmap::iterator pit = PM.begin();pit!=PM.end();pit++){
    gs_ret gr=getscoreandp( pit->first, tgp, args);
    gr.name=pit->first;
    tgscore=gr.gscore;
    if(args.needsort==false)
      gos<<pit->first<<"\t"<<PM[pit->first].size()<<"\t"<<tgscore<<"\t"<<gr.p_gsea<<"\t"<<gr.pval<<"\t"<<gr.pval<<"\t"<<0<<"\t"<<gr.gindex<<"\t"<<0<<endl;
    else{
      cerr<<"."<<flush;
      //vp.push_back(pair<double,pds> (tgp, pds(gr, pit->first) )); // a (pval,((score,index),name)) structure
      vp.push_back(gr);
    }
  }
  cerr<<endl;
  if(args.needsort){
    sort(vp.begin(),vp.end(),dpds_less_key());
    // calculate FDR
    vector<double> pvvec;
    for(int i=0;i<vp.size();i++){
      pvvec.push_back(vp[i].pval);
    }
    for(int i=0;i<pvvec.size();i++){
      pvvec[i]=pvvec[i]/(i+1)*pvvec.size();
    }
    if(pvvec[pvvec.size()-1]>1.0){
      pvvec[pvvec.size()-1]=1.0;
    }
    for(int i=pvvec.size()-2;i>=0;i--){
      if(pvvec[i]>pvvec[i+1]){
        pvvec[i]=pvvec[i+1];
      }
    }
    // output
    for(int i=0;i<vp.size();i++){
      // structure: (pval,((score,index),name))
      gs_ret &t_pr=vp[i];
      t_pr.fdr=pvvec[i];
      gos<<t_pr.name<<"\t"<<PM[t_pr.name].size()<<"\t"
	      <<t_pr.gscore<<"\t"<<t_pr.p_gsea<<"\t"<<t_pr.pval<<"\t"
	      <<t_pr.fdr<<"\t"<<i+1<<"\t"<<t_pr.gindex<<"\t"
	      <<0<<endl;
    }
  }//end if args.needsort
  if(args.outputfile!="-"){
    ofs.close();
  }
  return 0;
}

int main(int argc, char* argv[]){
  CallingArgs args;
  if(argc==1) {
    cout<<"mageckGSEA: A fast implementation of GSEA enrichment test. v"<<GSEA_VERSION<<endl;
    return 0;
  }
  parseArguments(argc,argv,args);
  
  init_vec(args); 
  
  double tgscore=0;
  double tgp=0.0;
  if(PM.count(args.pname)>0){
    gs_ret gr=getscoreandp( args.pname, tgp, args);
    tgscore=gr.gscore;

    //determine the type of the output
    ofstream ofs;
    if(args.outputfile!="-"){
      ofs.open(args.outputfile.c_str());
      if(!ofs.is_open()){  
        cerr<<"Error opening output file "<<args.outputfile<<endl;
        return -1;
      }
    }
    ostream &gos=(args.outputfile=="-")?cout:ofs;

    gos<<"Pathway\tSize\tES\tp\tp_permutation\tFDR\tRanking\tHits\tLFC\n";
    gos<<args.pname<<"\t"<<PM[args.pname].size()<<"\t"<<tgscore<<"\t"<<gr.p_gsea<<"\t"<<tgp<<"\t"<<tgp<<"\t"<<1<<"\t"<<gr.gindex<<"\t"<<0<<endl;
    
    if(args.outputfile!="-"){
      ofs.close();
    }

  }else{
    testallgseas(args);
  }
 
  end_vec(args);
  return 0;
}
