// Load the GraphGrind header files
// (the name gives away this evolved from Ligra)
#include "ligra-numa.h"

struct BFS_F
{
    // Variables (vertex properties)
    intT* Parents;

    // Auxiliary type definition
    struct cache_t {
        intT parents;
    };
    // Configuration
    static const bool use_cache = false;

    // Methods
    BFS_F(intT* _Parents) : Parents(_Parents) {}

    // Update method for single-thread usage
    inline bool update (intT s, intT d) {
        if(Parents[d] == -1)
        {
            Parents[d] = s;
            return 1;
        }
        else return 0;
    }
    // Thread-safe version of the update method
    inline bool updateAtomic (intT s, intT d) {
        return (CAS(&Parents[d],(intT)-1,s));
    }

    // Version of the update method helping the compiler to promote
    // intermediate values to registers
    inline void create_cache(cache_t &cache, intT d) {
        cache.parents=Parents[d];
    }
    inline bool update(cache_t &cache, intT s) {
        if (cache.parents==-1)
        {
            cache.parents=s;
            return 1;
        }
        else return 0;
    }
    inline void commit_cache(cache_t &cache, intT d) {
        Parents[d]=cache.parents;
    }

    // Convergence: cond function checks if vertex has been visited yet
    // Returns: true if vertex should still be processed; false if converged.
    inline bool cond (intT d) {
        return (Parents[d] == -1);
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex;
    const intT n = GA.n;
    const intT m = GA.m;
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();

    // Create vertex property: stores parent of each vertex. Initialized to -1,
    // except for start vertex.
    mmap_ptr<intT> Parents;
    Parents.part_allocate(part);
    map_vertexL(part, [&](intT j) {Parents[j]=-1;});
    Parents[start] = start;

    // Create the initial frontier consisting only of the start vertex
    partitioned_vertices Frontier = partitioned_vertices::create(
	n, start, GA.get_partition().V[start].getOutDegree());

    // Iterate until no more changes are made, signified by empty frontier
    int rounds = 0;
    while( !Frontier.isEmpty() ) {
        partitioned_vertices output
	    = edgeMap( GA, Frontier, BFS_F(Parents), m/20);
        Frontier.del();
        Frontier = output; //set new frontier
	++rounds;
    }
    Frontier.del();
    Parents.del();

    std::cerr << "Longest path: " << (rounds-1) << "\n";
}
