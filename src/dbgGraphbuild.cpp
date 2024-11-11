#include </home/zhangyicai/bcalm/gatb-core/gatb-core/src/gatb/debruijn/impl/GraphUnitigs.hpp>
#include "dbgGraphbuild.h"
#include <string>

/*void replaceSubstring(std::string& source, const std::string& target, const std::string& replacement) {
    size_t start_pos = source.find(target);
    if (start_pos != std::string::npos) {
        source.replace(start_pos, target.length(), replacement);
    }
}
bcalm_1::bcalm_1 ()  : Tool ("bcalm_1"){
    // old options, now using GATB's built-in options for kmer counting and graph creation
    // but TODO would be nice to integrate --nb-glue-partitions in gatb someday
/*	getParser()->push_back (new OptionOneParam ("-in", "input file",  true));
	getParser()->push_back (new OptionOneParam ("-out", "output prefix",  false, "unitigs"));
	getParser()->push_back (new OptionOneParam ("-k", "kmer size",  false,"31"));
	getParser()->push_back (new OptionOneParam ("-m", "minimizer size",  false,"8"));
	getParser()->push_back (new OptionOneParam ("-abundance", "abundance threshold",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-minimizer-type", "use lexicographical minimizers (0) or frequency based (1)",  false,"1"));
	getParser()->push_back (new OptionOneParam ("-dsk-memory", "max memory for kmer counting (MB)", false, "1500"));
	getParser()->push_back (new OptionOneParam ("-dsk-disk", "max disk space for kmer counting (MB)", false, "default"));

    // glue options
    getParser()->push_front (new OptionNoParam  ("--only-uf",   "(for debugging only) stop after UF construction", false));
	getParser()->push_front (new OptionNoParam  ("--uf-stats",   "display UF statistics", false));
	getParser()->push_back (new OptionOneParam ("--nb-glue-partitions", "number of glue files on disk",  false,"200"));


    IOptionsParser* graphParser = GraphUnitigsTemplate<32>::getOptionsParser(false);

    // hiding options
    if (IOptionsParser* p = graphParser->getParser(STR_KMER_ABUNDANCE_MIN_THRESHOLD))  {  p->setVisible(false); }
    if (IOptionsParser* p = graphParser->getParser(STR_HISTOGRAM_MAX))  {  p->setVisible(false); }
    if (IOptionsParser* p = graphParser->getParser(STR_SOLIDITY_KIND))  {  p->setVisible(false); } // oohh. multi-sample dbg construction someday maybe?
    if (IOptionsParser* p = graphParser->getParser(STR_URI_SOLID_KMERS))  {  p->setVisible(false); }
    
    // setting defaults
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_REPARTITION_TYPE)))  {  p->setDefaultValue ("1"); }
    if (Option* p = dynamic_cast<Option*> (graphParser->getParser(STR_MINIMIZER_TYPE)))  {  p->setDefaultValue ("1"); }

    getParser()->push_back(graphParser);


}
template <size_t span>
struct Functor  {  void operator ()  (bcalm_1 *bcalm)
    {
        typedef GraphUnitigsTemplate<span> GraphType;
        GraphType graph;


        graph = GraphType::create (bcalm->getInput(), false /* do not load unitigs after*/);

        // delete the .h5 file
        bool delete_h5_file = true;
        if (delete_h5_file)
        {
            // copies the h5 naming mechanism in GraphUnitigs.cpp
            string input = bcalm->getInput()->getStr(STR_URI_INPUT);
            string prefix;
            if (bcalm->getInput()->get(STR_URI_OUTPUT)) 
                prefix = bcalm->getInput()->getStr(STR_URI_OUTPUT);
            else
                prefix = System::file().getBaseName (input) + prefix;

            System::file().remove (prefix + ".h5");
        }

    }
};*/


void bcalm_1::execute (){

    std::cout << "BCALM 2, version " << 2; 
#ifdef GIT_SHA1
    std::cout << ", git commit " << GIT_SHA1;
#endif
    std::cout << std::endl;


    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** We launch Minia with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<Functor,bcalm_1*> (kmerSize, this);

}*/
void bcalmDbgGraphbuild(string illuminaFile,int kmerSize,int AbundanceMin,int cores){
    cout<<"--------------build graph-------------------------"<<endl;
	IBank *shortReads = Bank::open(illuminaFile);
    string solid_kmer_thr_str="1";
    string illuminaGraph = illuminaFile + "_k" + std::to_string(kmerSize) + ".h5";
    string illuminaFa = illuminaFile + "_k" + std::to_string(kmerSize) + "_a" + std::to_string(AbundanceMin) + ".fa";
	std::string arg="-kmer-size "+std::to_string(kmerSize)+" -out "+illuminaGraph+" -nb-cores "+std::to_string(cores)+" -abundance-min "+std::to_string(AbundanceMin);
    //int argc=5;
    cout<<"arg:"<<arg<<endl;
    cout<<"--------------end-------------------------"<<endl;
     /*char* argv[5]= {
        "./bcalm",
        "-in",
        illuminaFile.c_str(),
        "-kmer-size",
        std::to_string(kmerSize).c_str(),
        "-nb-cores",
        std::to_string(cores).c_str(),
        "",
        ""
    };
    
    if(argc > 1 && (   strcmp(argv[1],STR_VERSION)==0 || strcmp(argv[1],"-v")==0    )     ){
        std::cout << "BCALM 2, version " << 2;
     	std::cout << std::endl << "Using gatb-core version "<< System::info().getVersion() << std::endl;
	}

    try
    {
        // We run the tool with the provided command line arguments.
        bcalm_1().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
    }*/
    /*IProperties* props=bcalm_1 ().run(argc,argv);
    GraphType graph = GraphType::create (props, false);*/
    GraphType graph = GraphType::create (shortReads,arg.c_str());
	/*std::string prefix = illuminaFile;
    for (const auto &entry : filesystem::directory_iterator(".")) {
        if (entry.is_regular_file() && entry.path().filename().string().find(prefix) == 0) {
            std::string filename = entry.path().filename().string();
            filesystem::remove(entry.path());
		}
	}*/
	//filesystem::remove (illuminaGraph);
    //replaceSubstring(illuminaFile, ".fa", ".unitigs.fa");
    filesystem::copy("dummy.unitigs.fa",illuminaFa);
    //filesystem::remove (illuminaFile);
   

}
