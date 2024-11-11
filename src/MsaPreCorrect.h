#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <iostream>
#include <stack>
#include <map>

#include <malloc.h>
#include <getopt.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <gatb/gatb_core.hpp>
#include "mmpriv.h"
#include "minimap.h"
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "khash.h"



#define MIN_KMER_LEN 4
#define MIN_SOLID_THR 1

#define MIN_ERROR_RATE 0.0
#define MAX_ERROR_RATE 0.5

#define MIN_TRIALS 1
#define MAX_TRIALS 100

#define MIN_NB_BRANCH 1
#define MAX_NB_BRANCH 10000

#define MIN_THREADS 0
#define MAX_THREADS 64


// other constants
//#define MAX_READ_LEN 500000
#define MAX_PATH_LEN 1000

// default values for correction
#define DEF_ERROR_RATE 0.40
#define DEF_MAX_BRANCH 200
#define DEF_TRIALS 5
#define DEF_THREADS 0
#define DEF_MAX_ABUNDANCE 2147483647

// alignment scores
#define ALIGN_MATCH 1
#define ALIGN_MISMATCH -3
#define ALIGN_INDEL -2

// Constants for trimming and splitting
#define DEF_MIN_REGION_LG 100

// Constants for statistics
#define STAT_FOUND 0
#define STAT_FOUND_LEN1 1
#define STAT_TOOLONG 2
#define STAT_EXPLOSION 3
#define STAT_NOPATH 4

#define STAT_TAIL "TAIL"
#define STAT_END2END "END2END"
#define STAT_GAPEXTEND "GAPEXTEND"

// File extensions
#define HDF5_EXT ".h5"
#define FASTA_EXT ".fa"
#define FASTQ_EXT ".fq"

// Utilities 

typedef struct aln_ls_paf{ 
	uint32_t index;
	uint32_t mut;
	char qt;
	uint32_t rs;
	uint32_t re;
	uint32_t matchcount;
	uint32_t n_cigar;
	uint32_t *cigar;

}aln_ls_paf;

typedef struct {
	size_t n, m; 
	aln_ls_paf *a;
} aln_paf_v;

typedef struct {
	size_t n, m; 
	kstring_t *a;
} Short_reads;
void init_aln_ls_paf(aln_ls_paf* aln, uint32_t _index, char _qt, mm_reg1_t reg) {
	aln->index = _index;
	aln->rs=reg.rs;
	aln->re=reg.re;
    aln->qt = _qt;
    aln->n_cigar = reg.p->n_cigar;
    aln->cigar = (uint32_t *)malloc(sizeof(uint32_t) * aln->n_cigar);
    for (uint32_t i = 0; i < aln->n_cigar; ++i) {
        aln->cigar[i] = reg.p->cigar[i];
    }
}
// min macro
#define MIN(a,b) ((a) < (b) ? (a) : (b))

// max macro
#define MAX(a,b) ((a) > (b) ? (a) : (b))

// complement reverse a dna seq
void reverse(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    switch(source[len-i-1]) {
      case 'a':
        target[i] = 't';
        break;
      case 'c':
        target[i] = 'g';
        break;
      case 'g':
        target[i] = 'c';
        break;
      case 't':
        target[i] = 'a';
        break;
      case 'n':
        target[i] = 'n';
        break;
      case 'A':
        target[i] = 'T';
        break;
      case 'C':
        target[i] = 'G';
        break;
      case 'G':
        target[i] = 'C';
        break;
      case 'T':
        target[i] = 'A';
        break;
      case 'N':
        target[i] = 'N';
        break;
    }
  }
  target[len] = '\0';
}
//////////////////////////////////////////////////
// check wheter a file is present and readable 
bool is_readable( const std::string & file ) 
{ 
    std::ifstream fichier( file.c_str() ); 
    return !fichier.fail(); 
} 


typedef struct {
  int n,m;
  int* a;
}int_v;

typedef struct {
  int_v mapSid;
  std::map<std::string,int> i;
  std::map<std::string,int> d;
  std::map<std::string,int> x;
}ssnp;

typedef struct {
  int n,m;
  ssnp* snp=new ssnp[3];
}snp_v;


/*
判断是否有存储某个碱基，若存储过则将key键的value+1，若，没有存储过则insert一个键为key，value为1的键值对
*/
bool haveKey(std::map<std::string,int> &b,char* key){
  std::string base=key;
  auto search=b.find(base);
  if(search==b.end()){
    b.insert({base,1});
    return false; 
  }
  search->second=search->second+1;
  return true;  
}
bool find_mut(std::map<int,int> id_num,int sidx,int matchcount,float threshold){
  auto search=id_num.find(sidx);
  if(search!=id_num.end()){
    float num=search->second;
    float fre=num/matchcount;
    if(fre>threshold)
      return true;
  }
  return false;
}

bool add_map(int pos,int pos2,std::map<int,int_v> &map){
  auto search=map.find(pos);
  if(search==map.end()){
    int_v array;
    kv_init(array);
    kv_push(int,array,pos2);
    map.insert({pos,array});
    return false;
  } 
  int_v array=map[pos];
  kv_push(int,array,pos2);
  search->second=array;
  return true;
}

//记录短读sid存在snp的次数
bool add_mut(std::map<int,int> &mut,int sid){
  auto search=mut.find(sid);
  if(search==mut.end()){
    mut.insert({sid,1});
    return false;
  }
  search->second=search->second+1;
  return true;  
}

/*
根据每条short reads与long read比对的cigar字符串中的插入、删除、错配信息记录长读每个位置的snp
*/
void cigar_ite(int snp_len,char type,int &pos,int &srPos,std::map<int,ssnp> &pos_snp,char *lread,char* sread,int kmer_len){

    //M:匹配,I:插入,X:错配,D:删除      
  if(type=='=')//匹配且不存在变异，k-mer深度+1，并将长读和短读的指针后移snp_len个位置
  {
    pos=pos+snp_len;
    srPos=srPos+snp_len;
    return ;
  }

  //若没有return 表明该位置包含snp
  char* snp_base=new char[snp_len+1];
    //将变异碱基存在snp_base中
  for(int i=0;i<snp_len;i++){

    if(type=='I'||type=='X'){
      snp_base[i]=sread[srPos+i];
    }
    else{
      snp_base[i]=lread[srPos+i];
    }
  }  
  snp_base[snp_len]='\0';  
  /*若为插入则后移短读指针，删除则后移长读指针，错配同时移动短读和长读指针*/
  if(type=='I'){
    haveKey(pos_snp[pos].i,snp_base);
    srPos=srPos+snp_len;
  }
  else if(type=='D'){
    haveKey(pos_snp[pos].d,snp_base);
    pos=pos+snp_len;
  }
  else{
    haveKey(pos_snp[pos].x,snp_base);
    srPos=srPos+snp_len;
    pos=pos+snp_len;
  }
}   

/*
遍历b中的所有snp信息
*/
void find_max(std::map<string,int> &b,int &maxdepth,int &type,string &base,int &rsnp_len,int t,char* lread,int pos,int &num){
  if(b.empty()) return;
  
  for (auto it = b.begin(); it != b.end(); ++it) {
   
    string snp_base=it->first;
    int snp_len=snp_base.length();
    int depth=it->second;
    num+=depth;
    char* lbase=new char[snp_len+1];
    for(int i=0;i<snp_len;i++)//将pos位置的碱基存入lbase并和snp中的碱基比较
      lbase[i]=lread[pos+i];
    lbase[snp_len]='\0';  
    if(depth<5||depth<maxdepth) continue;
    else if((depth==maxdepth&&(type==2||type==0)&&snp_base.c_str()==lbase)||depth>maxdepth){//深度相等且变异碱基与长读相同或者深度更大
      maxdepth=depth;
      type=t;
      base=lbase;
      rsnp_len=snp_len;
    }
  }

}
void generate_mut(char* lread,int read_len,std::map<int,ssnp> pos_snp,std::map<int,int> &mutReads,std::map<int,int_v> map){
  /*
  遍历长读每个位置的snp信息
  */
  for(int pos=0;pos<read_len;pos++){
    int type=0;
    string base;
    int maxdepth=0,rsnp_len;
    int snp_num=0,map_num=0;
    int c=0;
    /*
    若pos位置存在变异，则遍历该位置存储的所有snp并存储深度最大的snp详细信息
    */
    if(pos_snp.find(pos)!=pos_snp.end()){
      find_max(pos_snp[pos].i,maxdepth,type,base,rsnp_len,0,lread,pos,snp_num);
      find_max(pos_snp[pos].d,maxdepth,type,base,rsnp_len,1,lread,pos,snp_num);
      find_max(pos_snp[pos].x,maxdepth,type,base,rsnp_len,2,lread,pos,snp_num);  

      /*
      计算在pos位置匹配上的短读数量
      */
      for(auto map_ite=map.begin();map_ite!=map.end();++map_ite){
        int start=map_ite->first;
        if(pos<start){break;}
        else{
          int_v end=map_ite->second;
          int i=0;
          for(int i=0;i<end.n;i++){
            if(pos<=end.a[i])
              map_num++;
          }
        }
      }
      c=map_num-snp_num;  
    }
    /*
    记录发生突变的reads
    */
    if(c>=5){
      int_v readsId=pos_snp[pos].mapSid;
      for (int i = 0; i <readsId.n; i++)
      {
        add_mut(mutReads,readsId.a[i]);
      } 
    }
  }
}
void update_lread(char* lread,int read_len,char* newLread,std::map<int,ssnp> pos_snp,std::map<int,int_v> map){
  int newpos=0;//新长读的初始指针
  /*
  遍历长读每个位置的snp信息
  */
  for(int pos=0;pos<read_len;pos++){
    int type=0;
    string base;
    int maxdepth=0,rsnp_len;
    int snp_num=0,map_num=0;
    int c=0;
    /*
    若pos位置存在变异，则遍历该位置存储的所有snp并存储深度最大的snp详细信息
    */
    if(pos_snp.find(pos)!=pos_snp.end()){
      find_max(pos_snp[pos].i,maxdepth,type,base,rsnp_len,0,lread,pos,snp_num);
      find_max(pos_snp[pos].d,maxdepth,type,base,rsnp_len,1,lread,pos,snp_num);
      find_max(pos_snp[pos].x,maxdepth,type,base,rsnp_len,2,lread,pos,snp_num);  
      for(auto map_ite=map.begin();map_ite!=map.end();++map_ite){
        int start=map_ite->first;
        if(pos<start){break;}
        else{
          int_v end=map_ite->second;
          int i=0;
          for(int i=0;i<end.n;i++){
            if(pos<=end.a[i])
              map_num++;
          }
        }
      }
      c=map_num-snp_num;  
    }
    if(c>=5||maxdepth<5){//若最大深度小于5或者该位置匹配上的短读数量大于5，不更新该位置的长读，直接将原read对应位置的碱基放入新read中
      newLread[newpos]=lread[pos];
      newpos++;
      continue;
    }

    if(type==0){
      for(int i=0;i<rsnp_len;i++){
        newLread[newpos]=base[i];
        newpos++;
      }
      newLread[newpos]=lread[pos];
      newpos++;
    }
    else if(type==1){
      pos=pos+rsnp_len-1;
    }
    else{
      for(int i=0;i<rsnp_len;i++){
        newLread[newpos]=base[i];
        newpos++;
      }
      pos=pos+rsnp_len-1;
    }
  }
  newLread[newpos]='\0';  
}
void update_snp(int lcigar,uint32_t* cigar,std::map<int,ssnp> &pos_snp, aln_ls_paf *ls_detail,int &pos,int &srPos,char *lread,char *sread,int kmer_len){
  for(int j=0;j<lcigar;++j){
    int snp_len=cigar[j]>>4;//cigar字符串的长度
    char type=MM_CIGAR_STR[cigar[j]&0xf];//变异的类型

    if(type!='='){//判断pos位置是否存在snp信息，没有则初始化一个
      if(pos_snp.find(pos)==pos_snp.end()){
        int_v v;
        kv_init(v);
        ssnp all_snp={v,{},{},{}};
        pos_snp[pos]=all_snp;
      }
      int_v mapid=pos_snp[pos].mapSid;
      kv_push(int,mapid,ls_detail->index);
      pos_snp[pos].mapSid=mapid;
    }
    //根据cigar中的信息更新pos_snp每个位置的snp信息，只需要snp信息
    cigar_ite(snp_len,type,pos,srPos,pos_snp,lread,sread,kmer_len);
  }
  
}
void prepare_snp(Sequence seq, int kmer_len,aln_paf_v alnPaf,std::map<int,int> &mutReads,mm_bseq1_t* sreads,char* newLread,float threshold){
  char *lread = seq.getDataBuffer();//长读序列
  int read_len = seq.getDataSize();
  int n_qReads=alnPaf.n;

  std::map<int,ssnp> pos_snp;
  std::map<int,ssnp> new_pos_snp;
  std::map<int,int_v> map;//所有长读匹配上的短读的起始位置map[0]={[0,249],[23,273],...}
  std::map<int,int_v> map1;
  int* mapNumber=new int[read_len];//在pos位置匹配上的短读数量

  /*分别处理每条短读的比对信息*/
  for(int i=0;i<n_qReads;i++){ 
    aln_ls_paf *ls_detail=&(alnPaf.a[i]);
    int lcigar=ls_detail->n_cigar;
    uint32_t* cigar=ls_detail->cigar;
    int srPos=0;//短读开始位置
    int pos=ls_detail->rs;//长读开始位置
    int pos2=ls_detail->re;
    char* sread=new char[sreads[ls_detail->index].l_seq];//短读序列
    
    auto search=map.find(pos);
    if(search==map.end()){
      int_v v;
      kv_init(v);
      kv_push(int,v,pos2);
      map.insert({pos,v});
    }
    else{
      int_v v=map[pos];
      kv_push(int,v,pos2);
      search->second=v;
    }
    if(ls_detail->qt=='-'){//考虑短读为反向互补链的情况
      char *buffer = sreads[ls_detail->index].seq;
      reverse(sread, buffer, sreads[ls_detail->index].l_seq);
    }
    else{
      sread=sreads[ls_detail->index].seq;
    }
    update_snp(lcigar,cigar,pos_snp,ls_detail,pos,srPos,lread,sread,kmer_len);
  }

  //根据pos_snp更新生成新的long reads
  generate_mut(lread,read_len,pos_snp,mutReads,map);
  int nomutAlin=0;
  for(int i=0;i<n_qReads;i++){ 
    aln_ls_paf* ls_detail=&(alnPaf.a[i]);
    bool is_mut = find_mut(mutReads, ls_detail->index, ls_detail->matchcount,threshold);
    if (is_mut) {
      nomutAlin++;
      continue;
    }
    
    int lcigar=ls_detail->n_cigar;
    uint32_t* cigar=ls_detail->cigar;
    int srPos=0;
    int pos=ls_detail->rs;
    int pos2=ls_detail->re;

    char* sread=new char[sreads[ls_detail->index].l_seq];

    auto search=map1.find(pos);
    if(search==map1.end()){
      int_v v;
      kv_init(v);
      kv_push(int,v,pos2);
      map1.insert({pos,v});
    }
    else{
      int_v v=map1[pos];
      kv_push(int,v,pos2);
      search->second=v;
    }

    if(ls_detail->qt=='-'){
      char *buffer = sreads[ls_detail->index].seq;
      reverse(sread, buffer, sreads[ls_detail->index].l_seq);
    }
    else{
      sread=sreads[ls_detail->index].seq;
    }
    update_snp(lcigar,cigar,new_pos_snp,ls_detail,pos,srPos,lread,sread,kmer_len);
  }
  cout<<"reads:"<<seq.getIndex()<<" rawAlnNum:"<<n_qReads<<" filterAlnNum:"<<nomutAlin<<endl;

  update_lread(lread,read_len,newLread,new_pos_snp,map1);

}
