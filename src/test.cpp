#include <iostream>
#include <gatb/gatb_core.hpp>
typedef GraphUnitigsTemplate<32> GraphType;


void bcalmDbgGraphbuild(string illuminaFile,int kmerSize,int AbundanceMin){
	IBank *shortReads = Bank::open(illuminaFile);
    string solid_kmer_thr_str="1";
    string illuminaGraph = illuminaFile + "_k" + std::to_string(kmerSize) + ".h5";
    string illuminaFa = illuminaFile + "_k" + std::to_string(kmerSize) + "_a" + std::to_string(AbundanceMin) + ".fa";
	std::string arg="-kmer-size "+std::to_string(kmerSize)+" -out "+illuminaGraph+" -abundance-min "+std::to_string(AbundanceMin);
    const char* argv = arg.c_str();
    GraphType graph = GraphType::create (shortReads,argv);
	std::string prefix = "dummy.unitigs.fa.glue";
    /*for (const auto &entry : filesystem::directory_iterator(".")) {
        if (entry.is_regular_file() && entry.path().filename().string().find(prefix) == 0) {
            std::string filename = entry.path().filename().string();
            std::filesystem::remove(entry.path());
		}
	}*/
	//System::file().remove (illuminaGraph);
    //System::file().copy("dummy.unitigs.fa",illuminaFa);
}
int main(){
    std::string illuminaFile="/home/zhangyicai/graphAligner/evalue/20strains/data/short_reads.fq";
    bcalmDbgGraphbuild(illuminaFile,21,1);
}