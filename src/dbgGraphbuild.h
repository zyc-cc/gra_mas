#ifndef dbgGraphbuild_h
#define dbgGraphbuild_h
#include <filesystem>
#include <gatb/gatb_core.hpp>
typedef GraphUnitigsTemplate<32> GraphType;

/*class bcalm_1 : public Tool
{
public:	
   
    // Constructor
    bcalm_1 ();

    // Actual job done by the tool is here
    void execute ();
    // Actual job done by the tool is here
};*/


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

#endif