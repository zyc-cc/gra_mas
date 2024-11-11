#ifndef Msa_h
#define Msa_h
#include <gatb/gatb_core.hpp>
#include "utils.hpp"
#include "MsaPreCorrect.h"
#include "map.h"


#define TO_STR2(x) #x
#define STR(x) TO_STR2(x)
#define LORDECSTRVERSION STR(LORDECVERSION)
#define GATBSTRVERSION STR(GATBVERSION)
//////////////////////////////////////////////////
// STACK class
class stack_element {
  public:

    Node node;
    int pos;
    char nt;

    stack_element(Node n, int p, char c);
};
stack_element::stack_element(Node n, int p, char c) {
  node=n;pos=p;nt=c;
}
// FLAGS
bool DEBUG = true;

// Global variables for correction parameters
float max_error_rate = DEF_ERROR_RATE;
int max_branch = DEF_MAX_BRANCH;
int max_trials = DEF_TRIALS;
int threads = DEF_THREADS;

int strict_mode = 0;

#ifndef GATB_V110
static const int KSIZE_4 = KMER_SPAN(3);
#endif


static const int progressBarWidth = 70;

std::string dirname(std::string source) {
  source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
  return source;
}

std::string basename(std::string pathname) {
  std::string dirn = dirname(pathname);
  size_t pos = pathname.find(dirn);
  size_t len = dirn.length();
  const std::string empty = "";
  return pathname.replace(pos, len, empty);
}

//////////////////////////////////////////////////

//////////////////////////////////////////////////
// string manipulation lowercase - uppercase
void copy_lower_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = tolower(source[i]);
  }
}

void copy_upper_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = toupper(source[i]);
  }
}
typedef struct {
	int n_processed, n_threads, n_fp;
	int64_t mini_batch_size;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;
	int n_parts;
	uint32_t *rid_shift;
	FILE *fp_split, **fp_parts;
  aln_paf_v* alnPafs;
  mm_bseq1_t *seq;
  int n_seq;
}pipeline_t1;

typedef struct {
	pipeline_t1 *p;
  int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg, *rep_len, *frag_gap;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
}step_t1;
static void *worker_pipeline1(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t1 *p = (pipeline_t1*)shared;
    if (step == 0) { // step 0: read sequences
      int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
      int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
      int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
      step_t1 *s;
      s = (step_t1*)calloc(1, sizeof(step_t1));
      if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
      else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
      //s->n_seq=p->n_seq;
      //s->seq=p->seq;
      if (s->seq) {
        s->p = p;
        for (i = 0; i < s->n_seq; ++i)
          s->seq[i].rid = p->n_processed++;
        s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
        for (i = 0; i < p->n_threads; ++i)
          s->buf[i] = mm_tbuf_init();
        s->n_reg = (int*)calloc(5 * s->n_seq, sizeof(int));
        s->seg_off = s->n_reg + s->n_seq; // seg_off, n_seg, rep_len and frag_gap are allocated together with n_reg
        s->n_seg = s->seg_off + s->n_seq;
        s->rep_len = s->n_seg + s->n_seq;
        s->frag_gap = s->rep_len + s->n_seq;
        s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
        for (i = 1, j = 0; i <= s->n_seq; ++i)
          if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
            s->n_seg[s->n_frag] = i - j;
            s->seg_off[s->n_frag++] = j;
            j = i;
          }
        return s;
		  } else free(s);
    } else if (step == 1) { // step 1: map
      if (p->n_parts > 0) merge_hits((step_t*)in);
      else kt_for(p->n_threads, worker_for, in, ((step_t1*)in)->n_frag);
      return in;
    } else if (step == 2) { // step 2: output
		  auto start = std::chrono::high_resolution_clock::now();
      void *km = 0;
      step_t1 *s = (step_t1*)in;
      aln_paf_v* alnPafs=p->alnPafs;
		  const mm_idx_t *mi = p->mi;
		  for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		  free(s->buf);
		  if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		  for (k = 0; k < s->n_frag; ++k) {
        int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
        for (i = seg_st; i < seg_en; ++i) {
          mm_bseq1_t *t = &s->seq[i];
          //cout<<"shortReads"<<i<<":"<<t->name<<endl;
          if (p->opt->split_prefix && p->n_parts == 0) { // then write to temporary files
            mm_err_fwrite(&s->n_reg[i],    sizeof(int), 1, p->fp_split);
            mm_err_fwrite(&s->rep_len[i],  sizeof(int), 1, p->fp_split);
            mm_err_fwrite(&s->frag_gap[i], sizeof(int), 1, p->fp_split);
            for (j = 0; j < s->n_reg[i]; ++j) {
              mm_reg1_t *r = &s->reg[i][j];
              mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
              if (p->opt->flag & MM_F_CIGAR) {
                mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
                mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
              }
            }
          } else if (s->n_reg[i] > 0) { // the query has at least one hit

            for (j = 0; j < s->n_reg[i]; ++j) {
              
              mm_reg1_t *r = &s->reg[i][j];
              
              assert(!r->sam_pri || r->id == r->parent);
              if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
                continue;
              int index = r->rid;
              //cout<<"longReads:"<<index<<endl;
              aln_ls_paf map_detail;
              char qt = "+-"[r->rev];
              init_aln_ls_paf(&map_detail, i, qt, s->reg[i][j]);
              std::mutex mtx;
              mtx.lock();
              kv_push(aln_ls_paf, alnPafs[index], map_detail);
              mtx.unlock();
            }
            //cout<<endl;
          } 
        }
        for (i = seg_st; i < seg_en; ++i) {
          for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
          free(s->reg[i]);
          free(s->seq[i].seq); free(s->seq[i].name);
          if (s->seq[i].qual) free(s->seq[i].qual);
          if (s->seq[i].comment) free(s->seq[i].comment);
        }
		  }
      free(s->reg); free(s->n_reg); free(s->seq); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
      km_destroy(km);
      if (mm_verbose >= 3)
        fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
      free(s);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> times = end - start;
      cout<<"print:"<<times.count()<<endl;
    }
    return 0;
}
int mm_map_file_frag1(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads, aln_paf_v* alnPafs,mm_bseq1_t* sReads,int n_seq)
{
	int i, pl_threads;
	pipeline_t1 pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t1));
	pl.n_fp = n_segs;
	pl.fp = open_bseqs(pl.n_fp, fn);
	if (pl.fp == 0) return -1;
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl.alnPafs=alnPafs;
  pl.seq=sReads;
  pl.n_seq=n_seq;
	if (opt->split_prefix)
		pl.fp_split = mm_split_init(opt->split_prefix, idx);
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline1, &pl, 3);
  
	free(pl.str.s);
	if (pl.fp_split) fclose(pl.fp_split);
	for (i = 0; i < pl.n_fp; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 1;
}

// target_file is long reads
// query_file is short reads
KSEQ_INIT(gzFile, gzread)
void runMinimap(const char *target_file,const char *query_file, aln_paf_v *alnPafs,std::map<int,int>* mutReads,int long_reads_num,mm_bseq1_t* sReads,int n_seq) {
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    int n_threads = 64;

    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(0, &iopt, &mopt);//初始化参数
    iopt.k=21;
    iopt.w=11;
    mopt.min_dp_max=60;
    mopt.min_chain_score = 30;
    mopt.min_cnt = 2;
    mopt.a =4, mopt.b =2;
    mopt.best_n = 2;
    mopt.end_bonus = 100;
    mopt.flag |= MM_F_ALL_CHAINS;
    mopt.flag |= MM_F_CIGAR; // perform alignment
    mopt.flag |= MM_F_EQX;
    mopt.flag |= MM_F_NO_LJOIN ;
    mopt.flag |= MM_F_SR  ;

    mm_idx_t *mi;

    kstring_t *s;
    void *km = 0;
    // // Open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(target_file, &iopt, 0);
    
    //*alnPafs = new aln_paf_v[long_reads_num];
    for (int i = 0; i < long_reads_num; ++i) {
      kv_init(alnPafs[i]);
    }
    auto start = std::chrono::high_resolution_clock::now();
    
    while ((mi = mm_idx_reader_read(r, n_threads)) != 0) {
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> times = end - start;
      cout<<"mm_idx_reader_read:"<<times.count()<<endl;
      mm_mapopt_update(&mopt, mi);

      mm_tbuf_t *tbuf = mm_tbuf_init();
      int short_reads_id=0;
      mm_map_file_frag1(mi, 1, &query_file, &mopt, n_threads, alnPafs,sReads,n_seq);
      //mm_map_file_frag(mi, 1, &query_file, &mopt, n_threads);
      end = std::chrono::high_resolution_clock::now();
      times = end - start;
      cout<<"mm_map_file_frag1:"<<times.count()<<endl;
      mm_tbuf_destroy(tbuf);
      mm_idx_destroy(mi);
      end = std::chrono::high_resolution_clock::now();
      times = end - start;
      cout<<"time:"<<times.count()<<endl;
      //fclose(file);
    }     
    
    mm_idx_reader_close(r); // 关闭索引读取器
    //kseq_destroy(ks);
    //gzclose(f);    
}




int msa(std::string pacbioFile,std::string illuminaFile,std::string outputFile,int thread,float threshold) {
  try {

    extern int opterr;
    opterr = 1;


    // parameters
    int
      kmer_len = MIN_KMER_LEN;
    std::string
      kmer_len_str = "4",
      solid_kmer_thr_str = "1",
      preCorrectedReadsfile,
      outTmpPath = "";


    //////// Count reads in input file
    BankFasta bsize(pacbioFile);
    BankFasta::Iterator itSeqSize(bsize);
    size_t max_read_len = 0;
    size_t seqSize;
    long long nbSeq = 0;
    for (itSeqSize.first(); !itSeqSize.isDone(); itSeqSize.next()) {
      seqSize = itSeqSize->getDataSize();
      if (seqSize > max_read_len) {
        max_read_len = seqSize;
      }
      nbSeq++;
    }
    max_read_len = max_read_len * 1.25;
    
    IBank *ptrBankPB = NULL;
    // Bank *ptrBankPB = NULL;
    try{
      // v106 simplified Bank interface (no BankRegistery anymore)// new in gatb-core-1.0.6
      ptrBankPB =  Bank::open(pacbioFile);
      
    }
    catch (gatb::core::system::Exception & bankPBExc){
      std::cout << "Error message PacBio bank " << bankPBExc.getMessage () << std::endl;
      return EXIT_FAILURE;
    }

    
    mm_bseq_file_t* fp= mm_bseq_open(illuminaFile.c_str());
    int short_reads_number;
    mm_bseq1_t* sReads=mm_bseq_read3(fp,500000000, 0,0,0, &short_reads_number);
   
    Iterator<Sequence> *itSeq= ptrBankPB->iterator();;

    // Create the output bank
    BankFasta output(outputFile); //argv[5]
    // Bank output(outReadFile); //argv[5]

    cout << "Found " << nbSeq << " reads.\n";
    cout << "Correcting reads...\n";

    ProgressManager *pmCorrect = new ProgressManager(nbSeq, "reads");
    std::map<int,int>* mutReads=new std::map<int,int>[nbSeq];
    aln_paf_v *alnPafs=new aln_paf_v[nbSeq];
    runMinimap(pacbioFile.c_str(),illuminaFile.c_str(),alnPafs,mutReads,nbSeq,sReads,short_reads_number);


    ISynchronizer *sync = System::thread().newSynchronizer();

    IDispatcher::Status status = Dispatcher(thread).iterate(itSeq, [&] (const Sequence& seq) {

      if (seq.getDataSize() >= kmer_len) {
        int read_len=seq.getDataSize();
        char *buffer = new char[max_read_len];
        char *newLread=new char[max_read_len];
        int *map_depth= NULL;

        if (seq.getDataSize() > max_read_len) {
          std::cout << "Too long read" << std::endl;
          exit(EXIT_FAILURE);
        }
        int long_idx=seq.getIndex();
        // Correct the read backward
        copy_upper_case(buffer, seq.getDataBuffer(), seq.getDataSize());
        buffer[seq.getDataSize()] = '\0';
        Sequence tmp_seq(buffer);
        tmp_seq._comment = seq.getComment();
        prepare_snp(tmp_seq, kmer_len,alnPafs[long_idx],mutReads[long_idx],sReads,newLread,threshold);
        Sequence s(newLread);
        s._comment = seq.getComment();
        {
          LocalSynchronizer local(sync);
          output.insert(s);
        }
        //cout<<"finish11"<<endl;
      
        delete [] buffer;
        delete [] newLread;
      }
    });
    delete [] alnPafs;

    cout<<"finish1"<<endl;

    std::cout << std::endl;
    delete pmCorrect;
    output.flush();
    delete ptrBankPB;
    return EXIT_SUCCESS;
  } catch (gatb::core::system::Exception & e){
    std::cout << "Error message " << e.getMessage () << std::endl;
    return EXIT_FAILURE;
  }
  
}


#endif