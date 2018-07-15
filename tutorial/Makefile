#define the long 64--bit for large scale graph
INTT = -DLONG 
CXXFLAGS = -fcilkplus -lcilkrts -O3 -DCILK=1 $(INTT) $(INTE) -I .. -I ../cilk-swan/include -L ../cilk-swan/lib -ldl  
HDR=papi_code.h utils.h IO.h parallel.h gettime.h quickSort.h parseCommandLine.h mm.h partitioner.h graph-numa.h ligra-numa.h
COMMON=$(patsubst %,../%,$(HDR))

ALL= BFS Components PageRank

#HILBERT=0, COO will use COO_CSR. For VEBO graph , COO_CSR is faster choice.
#LIBS_I_NEED= -DEDGES_HILBERT=1
LIBS_I_NEED= 

all: $(ALL)

#other option
CLIDOPT += -std=c++14
#CACHE collection, if PAPI_CACHE=1 collect and print the values
#CACHEOPT += -lpapi -DPAPI_CACHE=0
# not use atomic for forward and coo 
SEQOPT += -DPART96=1 
# NUMA special allocation - no NUMA on systems
#NUMAOPT += -DNUMA=1 -lnuma 
NUMAOPT += -DNUMA=0
#REUSE measurement, inital use reuse should not numa allocation
REUSEOPT += -DREUSE_DIST=0
#When using clang:
# OPT += -Wall -Wno-cilk-loop-control-var-modification

CXXFLAGS += $(NUMAOPT) $(CLIDOPT) $(SEQOPT) $(OPT) $(CACHEOPT) 

% : %.C $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS_I_NEED)

.PHONY : clean

clean :
	rm -f *.o $(ALL)
