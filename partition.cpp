// не протестированно 

#include "partition.h"
#include "random_chooser.h"

#include <unordered_set>
#include <algorithm>
#include <random>
#include <iostream>


// проверить, вдруг алгоритм leiden быстрее в libleidenalg TODO

static inline bool FlipCoin() {
    thread_local std::mt19937 engine(std::random_device{}());
    static std::bernoulli_distribution distribution(0.5);
    return distribution(engine);
}

size_t UniqueCnt(const igraph_vector_int_t* v) {
    std::unordered_set<int64_t> cnt;
    for (size_t i = 0; i < igraph_vector_int_size(v); ++i) {
        cnt.insert(v->stor_begin[i]);
    }
    return cnt.size();
}

Partition::Partition(const igraph_t* graph, double sample_fraction,
                    const igraph_vector_int_t* init_partition) 
  : n_nodes(igraph_vcount(graph))
  , n_edges(igraph_ecount(graph)) {
    this->_graph = graph; //destroy only in main()
    if (init_partition == NULL) {
        InitializePartition(sample_fraction);
    } else {
        igraph_vector_int_init_copy(&_membership, init_partition);
        n_comms = UniqueCnt(&_membership);
    }
    //std::cout << "partition " << n_nodes << " " << n_edges << "\n";
}

Partition::Partition() 
    : n_nodes(0)
    , n_edges(0)
    , n_comms(0)
    , _graph(NULL) {
        igraph_vector_int_init(&_membership, 0);
    }

Partition::Partition(const Partition& p) 
    : n_nodes(p.n_nodes)
    , n_edges(p.n_edges) {
        std::lock_guard<std::mutex> lock_other(p.mtx);
        this->_graph = p._graph;
        igraph_vector_int_init_copy(&_membership, &p._membership);
        n_comms = p.n_comms;
        _fittness = p._fittness;
        chooser = p.chooser;
}

Partition& Partition::operator=(const Partition& p) {
    if (this != &p) {
        std::scoped_lock lock(mtx, p.mtx);
        n_nodes = p.n_nodes;
        n_edges = p.n_edges;
        n_comms = p.n_comms;
        _fittness = p._fittness;
        igraph_vector_int_update(&_membership, &p._membership);
        _graph = p._graph;
        chooser = p.chooser;
    }
    return *this;
}

Partition::~Partition() {
    igraph_vector_int_destroy(&_membership);
}

void get_node_weights_for_modularity(const igraph_t* graph, igraph_vector_t* node_weights) {
    size_t n_nodes = igraph_vcount(graph);
    igraph_vector_int_t node_degrees; 
    igraph_vector_int_init(&node_degrees, n_nodes);
    igraph_degree(graph, &node_degrees, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_vector_init(node_weights, n_nodes);
    for (size_t i = 0; i < n_nodes; ++i) {
        node_weights->stor_begin[i] = node_degrees.stor_begin[i];
    }
    igraph_vector_int_destroy(&node_degrees);
}

void Partition::Optimize() {
    std::lock_guard<std::mutex> lock(mtx);
    igraph_vector_t node_weights; // сделать полем класса? TODO
    get_node_weights_for_modularity(_graph, &node_weights);
    if (igraph_vector_int_size(&_membership) != n_nodes) {
        igraph_vector_int_resize(&_membership, n_nodes);
    }
    igraph_community_leiden(_graph, NULL, &node_weights, 1.0/(2*n_edges), 
                    0.01, true, 3, &_membership, &n_comms, &_fittness);
    
    igraph_vector_destroy(&node_weights);
}

//Newman  O(|E|+|V|^2*steps)
//Walktrap O(|E||V|^2) in the worst case, O(|V|^2 log|V|) typically
//spinglass ?
//проблемка - в сишной версии нельзя указать, сколько кластеров я хочу, только максимальное количество итераций, а они не очень-то коррелируют
void Partition::NewmanSplit(igraph_vector_int_t* indices, igraph_integer_t comm_id_to_split) {
    igraph_t subgraph;
    igraph_induced_subgraph(_graph, &subgraph, igraph_vss_vector(indices), IGRAPH_SUBGRAPH_AUTO);
    igraph_vector_int_t* subgraph_membership = indices;
    igraph_community_leading_eigenvector(&subgraph, NULL, NULL, subgraph_membership, 5 /*steps ?*/, 
                                        NULL, NULL, false, NULL, NULL, NULL, NULL, NULL);
    for (size_t i = 0; i < igraph_vector_int_size(subgraph_membership); ++i) {
        subgraph_membership->stor_begin[i] = subgraph_membership->stor_begin[i] == 0 ? comm_id_to_split : n_comms;
    }
    
    for (size_t i = 0, j = 0; i < n_nodes; ++i) {
        if (_membership.stor_begin[i] == comm_id_to_split) {
            _membership.stor_begin[i] = subgraph_membership->stor_begin[j];
            j++;
        }
    }
    n_comms++;
    igraph_destroy(&subgraph);
}


void Partition::RandomSplit(igraph_vector_int_t* indices) {
    size_t size_to_split = std::max(size_t(1), size_t((std::rand() % igraph_vector_int_size(indices))/2)); // у Гала тут почему-то min, но выглядит как ошибка и нужно на самом деле max. но при этом я не уверена, что так вообще останутся связанные компоненты
    std::vector<size_t> indices_vector(igraph_vector_int_size(indices));
    for (size_t i = 0; i < indices_vector.size(); ++i) {
        indices_vector[i] = indices->stor_begin[i];
    }
    std::vector<size_t> indices_to_split = chooser.RandomChoice(indices_vector, size_to_split);
    for (const auto i : indices_to_split) {
        _membership.stor_begin[i] = n_comms;
    }
    n_comms++;
}

void copy_vector(igraph_vector_int_t* to, const std::vector<int32_t> &from) {
    igraph_vector_int_init(to, from.size());
    for (size_t i = 0; i < from.size(); ++i) {
        to->stor_begin[i] = from[i];
    }
}

void Partition::RandomMerge(size_t edges_subsample_size) {
    std::vector<int> candidate_edges_ids_vector = chooser.RandomChoice(n_edges, 
                                                         edges_subsample_size);
    igraph_vector_int_t candidate_edges_ids;
    copy_vector(&candidate_edges_ids, candidate_edges_ids_vector);
    igraph_eit_t candidate_edges;
    igraph_eit_create(_graph, igraph_ess_vector(&candidate_edges_ids), &candidate_edges);
    // не очень понятно, зачем по списку id делать итератор, чтобы из него потом снова получать id
    // надо проверить, возможно это обеспечивает проверку, что такое ребро есть, но в целом я это могу сделать
    // и сама, не аллоцируя лишнюю память TODO
    igraph_integer_t edge_id;
    igraph_integer_t from, to;
    igraph_integer_t comm1, comm2;
    while(!IGRAPH_EIT_END(candidate_edges)) {
        edge_id = IGRAPH_EIT_GET(candidate_edges);
        igraph_edge(_graph, edge_id, &from, &to);
        comm1 = _membership.stor_begin[from];
        comm2 = _membership.stor_begin[to];
        if (comm1 != comm2) {
            // if find an edge, connecting two different communities, 
            // merge these communities and break
            for (size_t i = 0; i < n_nodes; ++i) {
                if (_membership.stor_begin[i] == comm1) {
                    _membership.stor_begin[i] = comm2;
                } else if (_membership.stor_begin[i] == n_comms - 1) {
                    _membership.stor_begin[i] = comm1;
                }
            }
            n_comms--;
            break;
        }
        IGRAPH_EIT_NEXT(candidate_edges);
    }
    igraph_eit_destroy(&candidate_edges);
    igraph_vector_int_destroy(&candidate_edges_ids);
}

void Partition::Mutate() {
    std::lock_guard<std::mutex> lock(mtx);
    if (FlipCoin()) {
        // split a random community
        igraph_integer_t comm_id_to_split = std::rand() % n_comms;
        std::vector<int> indices_to_split_vector;
        for (size_t i = 0; i < n_nodes; ++i) {
            if (_membership.stor_begin[i] == comm_id_to_split) {
                indices_to_split_vector.push_back(i);
            }
        }
        igraph_vector_int_t indices_to_split;
        copy_vector(&indices_to_split, indices_to_split_vector);
        if (igraph_vector_int_size(&indices_to_split) > 2) {
            size_t min_comm_size_newman = 10; // hardcoded value
            if (igraph_vector_int_size(&indices_to_split) > min_comm_size_newman) {
                if (FlipCoin()) {
                    // split via Newman's leading eigenvector method for detecting community structure
                    NewmanSplit(&indices_to_split, comm_id_to_split);
                } else {
                    // split randomly 
                    RandomSplit(&indices_to_split);
                }
            } else {
                // split randomly
                RandomSplit(&indices_to_split);
            }
        } 
        igraph_vector_int_destroy(&indices_to_split);

    } else {
        // randomly merge two connected communities 
        size_t edges_subsample_size = 10; // hardcoded value
        RandomMerge(edges_subsample_size);
    }
}

void GetVertexIdVector(const igraph_t* graph, igraph_vector_int_t* out) {
    igraph_vit_t vit;
    igraph_vit_create(graph, igraph_vss_all(), &vit);
    igraph_vector_int_init(out, igraph_vcount(graph));
    size_t i = 0;
    while (!IGRAPH_VIT_END(vit)) {
        out->stor_begin[i] = IGRAPH_VIT_GET(vit);
        IGRAPH_VIT_NEXT(vit);
        i++;
    }
    igraph_vit_destroy(&vit);
}

void Partition::InitializePartition(double sample_fraction) {
    igraph_t subgraph;
    if (FlipCoin()) {
        // sample nodes
        auto subsample = chooser.RandomChoice(n_nodes, sample_fraction*n_nodes);
        igraph_vector_int_t subsample_nodes_ids; 
        copy_vector(&subsample_nodes_ids, subsample);
        igraph_induced_subgraph(_graph, &subgraph, igraph_vss_vector(&subsample_nodes_ids), 
                                                    IGRAPH_SUBGRAPH_AUTO);
        igraph_vector_int_destroy(&subsample_nodes_ids);
    } else {
        // sample edges
        auto subsample = chooser.RandomChoice(n_edges, sample_fraction*n_edges);
        igraph_vector_int_t subsample_edges_ids; 
        copy_vector(&subsample_edges_ids, subsample);
        igraph_subgraph_from_edges(_graph, &subgraph, igraph_ess_vector(&subsample_edges_ids), 
                                                                                    true);
        igraph_vector_int_destroy(&subsample_edges_ids);
    }

    igraph_vector_t node_weights;
    get_node_weights_for_modularity(&subgraph, &node_weights);
    auto n_nodes_subgraph = igraph_vcount(&subgraph);
    igraph_vector_int_t subgraph_membership; 
    igraph_vector_int_init(&subgraph_membership, n_nodes_subgraph);
    igraph_community_leiden(&subgraph, NULL, &node_weights, 1.0/(2*igraph_ecount(&subgraph)), 
            0.01, false, 2, &subgraph_membership, &n_comms, &_fittness); // objective function - modularity
    
    igraph_vector_destroy(&node_weights);

    igraph_vector_int_t subgraph_vertices;
    GetVertexIdVector(&subgraph, &subgraph_vertices);
    igraph_vector_int_init(&_membership, n_nodes);
    for (size_t i = 0; i < n_nodes; ++i) {
        _membership.stor_begin[i] = -1;
    }
    
    for (size_t i = 0; i < n_nodes_subgraph; ++i) {
        _membership.stor_begin[subgraph_vertices.stor_begin[i]] = subgraph_membership.stor_begin[i];
    }

    igraph_vector_int_destroy(&subgraph_vertices);
    igraph_vector_int_destroy(&subgraph_membership);
    igraph_destroy(&subgraph);

    for (size_t i = 0; i < n_nodes; ++i) {
        if (_membership.stor_begin[i] == -1) {
            _membership.stor_begin[i] = n_comms++;
        }
    }
}
