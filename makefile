PLATFORM=$(shell uname -s)
GPP=$(CXX)
GATB=/home/zhangyicai/bcalm/gatb-core/gatb-core
GATB_LIB=$(GATB)/build/lib
CPPFLAGS=-Wall -Wextra -std=c++17 -O3 -Iconcurrentqueue -Izstr/src -Iparallel-hashmap/parallel_hashmap/ `pkg-config --cflags protobuf` `pkg-config --cflags libsparsehash` -Wno-unused-parameter -IMEMfinder/src -I`jemalloc-config --includedir` -Iminimap2/ -I$(GATB)/src  -I$(GATB)/build/include -I$(GATB)/thirdparty -L$(GATB_LIB)   -Lminimap2/ -lminimap2 -lgatbcore -lhdf5 -ldl -lz -lpthread

ODIR=obj
BINDIR=bin
SRCDIR=src

LIBS=-lm -lz -lboost_program_options `pkg-config --libs protobuf` -lsdsl
JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`

_DEPS =vg.pb.h fastqloader.h dbgGraphbuild.h Msa.h MsaPreCorrect.h utils.hpp GraphAlignerWrapper.h BigraphToDigraph.h stream.hpp Aligner.h ThreadReadAssertion.h AlignmentGraph.h CommonUtils.h GfaGraph.h ReadCorrection.h MinimizerSeeder.h AlignmentSelection.h EValue.h MEMSeeder.h DNAString.h DiploidHeuristic.h 
DEPS = $(patsubst %, $(SRCDIR)/%, $(_DEPS))

_OBJ =utils.o Aligner.o vg.pb.o fastqloader.o BigraphToDigraph.o ThreadReadAssertion.o AlignmentGraph.o CommonUtils.o GraphAlignerWrapper.o GfaGraph.o ReadCorrection.o MinimizerSeeder.o AlignmentSelection.o EValue.o MEMSeeder.o DNAString.o DiploidHeuristic.o
OBJ = $(patsubst %, $(ODIR)/%, $(_OBJ))

ifeq ($(PLATFORM),Linux)
   JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -Wl,-Bstatic -ljemalloc -Wl,-Bdynamic `jemalloc-config --libs`
   LINKFLAGS = $(CPPFLAGS) -Wl,-Bstatic $(LIBS) -Wl,-Bdynamic -Wl,--as-needed -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS) `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`  -I$(GATB)/include  -L$(GATB)/lib -lgatbcore -lhdf5 -ldl -lz -lpthread
else
   CPPFLAGS += -D_LIBCPP_DISABLE_AVAILABILITY
   JEMALLOCFLAGS= -L`jemalloc-config --libdir` -Wl,-rpath,`jemalloc-config --libdir` -ljemalloc `jemalloc-config --libs`
   LINKFLAGS = $(CPPFLAGS) $(LIBS) -lpthread -pthread -static-libstdc++ $(JEMALLOCFLAGS) `pkg-config --libs libdivsufsort` `pkg-config --libs libdivsufsort64`
endif

VERSION := Branch $(shell git rev-parse --abbrev-ref HEAD) commit $(shell git rev-parse HEAD) $(shell git show -s --format=%ci)

$(shell mkdir -p bin)
$(shell mkdir -p obj)

$(BINDIR)/GraphAligner: $(ODIR)/AlignerMain.o $(OBJ) MEMfinder/lib/memfinder.a
	$(GPP) -g -o $@ $^ $(LINKFLAGS) 

$(ODIR)/GraphAlignerWrapper.o: $(SRCDIR)/AlignerMain.cpp  $(SRCDIR)/GraphAligner.h $(SRCDIR)/NodeSlice.h $(SRCDIR)/WordSlice.h $(SRCDIR)/ArrayPriorityQueue.h $(SRCDIR)/ComponentPriorityQueue.h $(SRCDIR)/GraphAlignerVGAlignment.h $(SRCDIR)/GraphAlignerGAFAlignment.h $(SRCDIR)/GraphAlignerBitvectorBanded.h $(SRCDIR)/GraphAlignerBitvectorCommon.h $(SRCDIR)/GraphAlignerCommon.h  $(DEPS)

$(ODIR)/AlignerMain.o: $(SRCDIR)/AlignerMain.cpp  $(DEPS) $(OBJ)
	$(GPP) -g -c -o $@ $< $(CPPFLAGS) -DVERSION="\"$(VERSION)\""

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(GPP) -g -c -o $@ $< $(CPPFLAGS)

$(ODIR)/vg.pb.o: $(SRCDIR)/vg.pb.cc
	$(GPP) -g -c -o $@ $< $(CPPFLAGS)

$(SRCDIR)/%.pb.cc $(SRCDIR)/%.pb.h: $(SRCDIR)/%.proto
	protoc -I=$(SRCDIR) --cpp_out=$(SRCDIR) $<

MEMfinder/lib/memfinder.a:
	$(MAKE) -C MEMfinder lib DEBUGFLAG="-DNDEBUG"

all: $(BINDIR)/GraphAligner

clean:
	rm -f $(ODIR)/*
	rm -f $(BINDIR)/*
	rm -f $(SRCDIR)/vg.pb.cc
	rm -f $(SRCDIR)/vg.pb.h
	$(MAKE) -C MEMfinder clean