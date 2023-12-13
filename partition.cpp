#include "partition.h"
#include "random_chooser.h"

#include <set>

//igraph_vector_reserve(igraph_vector_t* v, igraph_integer_t capacity);

bool FlipCoin() {
    return static_cast<double>(std::rand()) / RAND_MAX > 0.5;
}

int UniqueCnt(const igraph_vector_int_t* v) {
    set<int> cnt;
    for (size_t i = 0; i < igraph_vector_int_size(v); ++i) {
        cnt.insert(VECTOR(v)[i]);
    }
    return cnt.size();
}

Partition::Partition(igraph_t* graph, double sample_fraction,
                    igraph_vector_int_t* init_partition) 
    : n_nodes(igraph_vcount(graph))
    , n_edges(igraph_ecount(graph)) 
    , sample_size_nodes(n_nodes * sample_fraction)
    , sample_size_edges(n_edges * sample_fraction) {
        this->_graph = graph; //возможно стоит всё-таки делать поверхностную копию и destroy тлько в main
        if (init_partition == NULL) {
            InitializePartition();
        } else {
            igraph_vector_int_update(_membership, init_partition);
            n_comms = UniqueCnt(_membership);
        }
    }

Partition::~Partition() {
    igraph_destroy((igraph_t*)this->_graph);
    delete this->_graph;
    igraph_vector_int_destroy(_membership);
    delete this->_membership;
}

Partition Partition::Optimize() {

}

void Partition::NewmanSplit(const vector<size_t>& indices, 
    size_t /*не уверена насчёт типа данных*/ comm_id_to_split) {

    }


void Partition::RandomSplit(const vector<size_t>& indices) {

}

Partition Partition::Mutate() {

}

// есть некая вероятность, что нужен только igraph, потому что в igraph/src/community/leiden.c есть igraph_community_leiden,
// и Гал использует именно его. Но, возможно, там медленнее, надо проверять
void Partition::InitializePartition() {
    igraph_t subgraph;
    if (FlipCoin()) {
        // sample nodes
        auto subsample = chooser.RandomChoice(n_nodes, sample_size_nodes);
        igraph_vector_int_t subsample_nodes_ids;
        igraph_vector_int_view(&subsample_nodes_ids, subsample.data(), subsample.size());
        igraph_induced_subgraph(_graph, &subgraph, igraph_vss_vector(&subsample_nodes_ids), 
                                                    IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH);
    } else {
        // sample edges
        auto subsample = chooser.RandomChoice(n_edges, sample_size_edges);
        igraph_vector_int_t subsample_edges_ids;
        igraph_vector_int_view(&subsample_edges_ids, subsample.data(), subsample.size());
        igraph_subgraph_from_edges(_graph, &subgraph, igraph_ess_vector(&subsample_edges_ids), 
                                                                                    true);
    }

    // Leiden algorithm on subgraph
    // https://igraph.org/c/doc/igraph-Community.html#igraph_community_leiden
    // igraph/examples/simple/igraph_community_leiden.c 

    igraph_vector_int_t node_weights;
    igraph_vector_init(&node_weights, igraph_vcount(&subgraph));
    igraph_degree(&subgraph, &node_weights, igraph_vss_all(), IGRAPH_ALL, true);

    igraph_integer_t nb_clusters;
    igraph_real_t quality;
    igraph_community_leiden(&subgraph, NULL, &node_weights, 1.0/(2*igraph_ecount(&subgraph)), 
                            0.01, false, 2, _membership, &nb_clusters, &quality); //это настройки для того, чтобы objective function была modularity
    


    igraph_destroy(&subgraph);
}