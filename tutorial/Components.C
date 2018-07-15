// Load the GraphGrind header files
// (the name gives away this evolved from Ligra)
#include "ligra-numa.h"

struct CC_F
{
    // Variables (vertex properties)
    intT* IDs;
    intT* prevIDs;

    // Auxiliary type definition
    struct cache_t { };
    // Configuration
    static const bool use_cache = false; // disable the create/update/commit cache methods

    // Methods
    CC_F(intT* _IDs, intT* _prevIDs) :
        IDs(_IDs), prevIDs(_prevIDs) {}

    // Update method for single-thread usage
    inline bool update(intT s, intT d) {
	return updateAtomic( s, d ); // for your convenience
    }
    // Thread-safe version of the update method
    inline bool updateAtomic (intT s, intT d) {
	// TODO ...
    }

    // Version of the update method helping the compiler to promote
    // intermediate values to registers
    inline void create_cache(cache_t &cache, intT d) { }
    inline bool update(cache_t &cache, intT s) { return 0; }
    inline void commit_cache(cache_t &cache, intT d) { }

    // Convergence: cond function checks if vertex has been visited yet
    // Returns: true if vertex should still be processed; false if converged.
    inline bool cond (intT d) {
	// TODO ...
    }
};

// function used by VertexMap to sync prevIDs with IDs
struct CC_Vertex_F
{
    // Variables (vertex properties)
    intT* IDs;
    intT* prevIDs;

    // Constructor of the class
    CC_Vertex_F(intT* _IDs, intT* _prevIDs) :
        IDs(_IDs), prevIDs(_prevIDs) {}

    // Called for every vertex i recorded in the frontier
    inline bool operator () (intT i) {
        prevIDs[i] = IDs[i];
        return 1;
    }
};

// To construct a summary histogram of component sizes
struct CC_Vertex_Count
{
    const intT* IDs;
    CC_Vertex_Count( const intT* _IDs ) : IDs(_IDs) {}
    inline bool operator () ( intT i ) {
        return IDs[i] == i;
    }
};

template <class GraphType>
void Compute(GraphType &GA, long start)
{
    typedef typename GraphType::vertex_type vertex;
    const partitioner &part = GA.get_partitioner();
    const int perNode = part.get_num_per_node_partitions();
    intT n = GA.n;
    intT m = GA.m;

    // Create vertex properties: stores current and previous IDs.
    // Initialized to vertex ID.
    mmap_ptr<intT> IDs;
    IDs.part_allocate (part);
    mmap_ptr<intT> prevIDs;
    prevIDs.part_allocate (part);
    map_vertexL(part,[&](intT j) {IDs[j]=j;});

    // Create the initial frontier, consisting of all vertices
    partitioned_vertices Frontier = partitioned_vertices::bits(part,n,m);

    // Iterate until no more changes are made, signified by empty frontier
    int rounds = 0;
    while(!Frontier.isEmpty()) {
        vertexMap(part,Frontier,CC_Vertex_F(IDs,prevIDs));
        partitioned_vertices output
	    = edgeMap(GA, Frontier, CC_F(IDs,prevIDs), m/20);
        Frontier.del();
        Frontier = output;
	++rounds;
    }

    // Now calculate the number of components.
    intT ncomponents = 0;
    map_vertexL( part, [&]( intT v ) {
	    if( IDs[v] == v )
		__sync_fetch_and_add( &ncomponents, 1 );
	} );

    std::cerr << "Required " << rounds << " rounds; found "
	      << ncomponents << " components\n";

    Frontier.del();
    IDs.del();
    prevIDs.del();
}

