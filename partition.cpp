// не протестированно 

#include "partition.h"
#include "random_chooser.h"

#include <set>

// всякие перкладывания значений из одного вектора в другой точно надо параллелить TODO
// хотя и возникнет много вложенной параллелизации, которой нет

// проверить, вдруг алгоритм leiden быстрее в libleidenalg TODO

// не использовать VECTOR, а просто тыкать stor_begin, наверное это быстрее TODO
// то же самое касается всех макросов для итераторов
// for (double *ptr = vec.stor_begin; ptr != vec.end; ++ptr) { ... } - но авторы igraph не ракомендуют так делать, мол библиотека и поменяться может

bool FlipCoin() {
    return static_cast<double>(std::rand()) / RAND_MAX > 0.5;
}

size_t UniqueCnt(const igraph_vector_int_t* v) {
    set<int> cnt; // не int, а что там у них
    for (size_t i = 0; i < igraph_vector_int_size(v); ++i) {
        cnt.insert(VECTOR(v)[i]);
    }
    return cnt.size();
}

Partition::Partition(igraph_t* graph, double sample_fraction,
                    const igraph_vector_int_t* init_partition) 
  : n_nodes(igraph_vcount(graph))
  , n_edges(igraph_ecount(graph)) {
    this->_graph = graph; //destroy only in main()
    if (init_partition == NULL) {
        InitializePartition(sample_fraction);
    } else {
        igraph_vector_int_update(&_membership, init_partition);
        n_comms = UniqueCnt(&_membership);
    }
}

Partition::~Partition() {
    delete this->_graph;
    igraph_vector_int_destroy(&_membership);
}

Partition Partition::Optimize() {
    igraph_vector_int_t node_weights; // сделать полем класса, 
    // чтобы не аллоцировать память каждый раз? пока не очень понимаю, 
    // сколько раз вызывается для 1 объекта. если больше 1, то имеет смысл TODO
    igraph_vector_int_init(&node_weights, n_nodes);
    igraph_degree(_graph, &node_weights, igraph_vss_all(), IGRAPH_ALL, true);

    igraph_community_leiden(_graph, NULL, &node_weights, 1.0/(2*n_edges), 
                    0.01, true, 3, &_membership, &n_comms, &_fittness);
    
    igraph_vector_int_destroy(&node_weights);
    return this;
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
        VECTOR(subgraph_membership)[i] = VECTOR(*subgraph_membership)[i] == 0 ? comm_id_to_split : n_comms;
    }
    
    for (size_t i = 0, j = 0; i < n_nodes; ++i) {
        if (VECTOR(_membership)[i] == comm_id_to_split) {
            VECTOR(_membership)[i] = VECTOR(*subgraph_membership)[j];
            j++;
        }
    }
    n_comms++;
    igraph_destroy(&subgraph);
}


void Partition::RandomSplit(igraph_vector_int_t* indices) {
    size_t size_to_split = max(1, (std::rand() % igraph_vector_int_size(indices))/2); // у Гала тут почему-то min, но выглядит как ошибка
    std::vector<size_t> indices_vector(igraph_vector_int_size(indices));
    for (size_t i = 0; i < indices_vector.size(); ++i) {
        indices_vector[i] = VECTOR(*indices)[i];
    }
    std::vector<size_t> indices_to_split = chooser.RandomChoice(indices_vector, size_to_split);
    for (const auto i : indices_to_split) {
        VECTOR(_membership)[i] = n_comms;
    }
    /*это ошибка, 
    и нужно на самом деле max. но при этом я не уверена, что так вообще 
    останутся связанные компоненты
    */
    n_comms++;
}


void Partition::RandomMerge(size_t edges_subsample_size) {
    std::vector<int> candidate_edges_ids_vector = chooser.RandomChoice(n_edges, 
                                                         edges_subsample_size);
    igraph_vector_int_t candidate_edges_ids;
    igraph_vector_int_view(&candidate_edges_ids, candidate_edges_ids_vector.data(),
                                candidate_edges_ids_vector.size());
    igraph_eit_t candidate_edges;
    igraph_eit_create(_graph, igraph_ess_vector(candidate_edges_ids), &candidate_edges);
    // не очень понятно, зачем по списку id делать итератор, чтобы из него потом снова получать id
    // надо проверить, возможно это обеспечивает проверку, что такое ребро есть, но в целом я это могу сделать
    // и сама, не аллоцируя лишнюю память TODO
    igraph_integer_t edge_id;
    igraph_integer_t from, to;
    igraph_integer_t comm1, comm2;
    while(!IGRAPH_EIT_END(candidate_edges)) {
        edge_id = IGRAPH_EIT_GET(candidate_edges);
        igraph_edge(_graph, edge_id, &from, &to);
        comm1 = VECTOR(_membership)[from];
        comm2 = VECTOR(_membership)[to];
        if (comm1 != comm2) {
            // if find an edge, connecting two different communities, 
            // merge these communities and break
            for (size_t i = 0; i < n_nodes; ++i) {
                if (VECTOR(_membership)[i] == comm1) {
                    VECTOR(_membership)[i] = comm2;
                } else if (VECTOR(_membership)[i] == n_comms - 1) {
                    VECTOR(_membership)[i] = comm1;
                }
            }
            n_comms--;
            break;
        }
        IGRAPH_EIT_NEXT(candidate_edges);
    }
    igraph_eit_destroy(&candidate_edges);
}

Partition Partition::Mutate() {
    if (FlipCoin()) {
        // split a random community
        igraph_integer_t comm_id_to_split = std::rand() % n_comms;
        std::vector<size_t> indices_to_split_vector;
        for (size_t i = 0; i < n_nodes; ++i) {
            if (VECTOR(*_membership)[i] == comm_id_to_split) {
                indices_to_split_vector.push_back(i);
            }
        }
        igraph_vector_int_t indices_to_split;
        igraph_vector_int_view(indices_to_split, indices_to_split_vector.data(), 
                                                indices_to_split_vector.size());
        if (igraph_vector_int_size(indices_to_split) > 2) {
            size_t min_comm_size_newman = 10; // hardcoded value
            if (igraph_vector_int_size(indices_to_split) > min_comm_size_newman) {
                if (FlipCoin()) {
                    // split via Newman's leading eigenvector method for detecting community structure
                    NewmanSplit(indices_to_split, comm_id_to_split);
                } else {
                    // split randomly 
                    RandomSplit(indices_to_split);
                }
            } else {
                // split randomly
                RandomSplit(indices_to_split);
            }
        } 

    } else {
        // randomly merge two connected communities 
        size_t edges_subsample_size = 10; // hardcoded value
        RandomMerge(edges_subsample_size);
    }
    return this;
}


//перенести функцию внутрб igraph либо переименовать под свой нейминг TODO
void igraph_vector_int_init_with(igraph_vector_int_t* v, size_t size, int fill) {
    igraph_vector_int_init(v, size);
    for (size_t i = 0; i < size; ++i) {
        VECTOR(v)[i] = fill;
    }
}

void GetVertexIdVector(const igraph_t* graph, igraph_vector_int_t* out) {
    igraph_vit_t vit;
    igraph_vit_create(graph, igraph_vss_all(), &vit);
    igraph_vector_int_init(out, igraph_vcount(graph));
    size_t i = 0;
    while (!IGRAPH_VIT_END(vit)) {
        VECTOR(out)[i] = IGRAPH_VIT_GET(vit);
        IGRAPH_VIT_NEXT(vit);
        i++;
    }
    igraph_vit_destroy(vit);
}

void Partition::InitializePartition(double sample_fraction) {
    igraph_t subgraph;
    if (FlipCoin()) {
        // sample nodes
        auto subsample = chooser.RandomChoice(n_nodes, sample_fraction*n_nodes);
        igraph_vector_int_t subsample_nodes_ids; //это вьюшка, destroy не нужен
        igraph_vector_int_view(&subsample_nodes_ids, subsample.data(), subsample.size());
        igraph_induced_subgraph(_graph, &subgraph, igraph_vss_vector(&subsample_nodes_ids), 
                                                    IGRAPH_SUBGRAPH_AUTO);
    } else {
        // sample edges
        auto subsample = chooser.RandomChoice(n_edges, sample_fraction*n_edges);
        igraph_vector_int_t subsample_edges_ids; //это вьюшка, destroy не нужен
        igraph_vector_int_view(&subsample_edges_ids, subsample.data(), subsample.size());
        igraph_subgraph_from_edges(_graph, &subgraph, igraph_ess_vector(&subsample_edges_ids), 
                                                                                    true);
    }

    auto n_nodes_subgraph = igraph_vcount(&subgraph);
    igraph_vector_int_t node_weights;
    igraph_vector_int_init(&node_weights, n_nodes_subgraph);
    igraph_degree(&subgraph, &node_weights, igraph_vss_all(), IGRAPH_ALL, true);

    igraph_vector_int_t subgraph_membership; 
    igraph_vector_int_init(&subgraph_membership, n_nodes_subgraph);
    igraph_community_leiden(&subgraph, NULL, &node_weights, 1.0/(2*igraph_ecount(&subgraph)), 
            0.01, false, 2, &subgraph_membership, &n_comms, &_fittness); // objective function - modularity
    
    igraph_vector_int_destroy(&node_weights);

    igraph_vector_int_t subgraph_vertices;
    GetVertexIdVector(&subgraph, subgraph_vertices);
    igraph_vector_int_init_with(&_membership, n_nodes, -1);
    
    for (size_t i = 0; i < n_nodes_subgraph; ++i) {
        VECTOR(_membership)[VECTOR(subgraph_vertices)[i]] = VECTOR(subgraph_membership)[i];
    }

    igraph_vector_int_destroy(&subgraph_vertices);
    igraph_vector_int_destroy(&subgraph_membership);
    igraph_destroy(&subgraph);

    for (size_t i = 0; i < n_nodes; ++i) {
        if (VECTOR(_membership)[i] == -1) {
            VECTOR(_membership)[i] = n_comms++;
        }
    }
}
