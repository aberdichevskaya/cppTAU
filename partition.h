#pragma once

#include "random_chooser"
#include "igraph/include/igraph.h"

#include <ctime>

#include <vector>

// https://igraph.org/c/doc 


// мб стоит брать не igraph, а Graph из libleidenalg/include/GraphHelper.h
// с другой стороны, subgraph там кажется нет

// может быть написать graph_helper как в libleidenalg?
// не пойму, надо ли хранить указатель на igraph, или норм хранить объект

class Partition {
 public:
    Partition(const igraph_t* G_ig, double sample_fraction=0.5, 
                igraph_vector_int_t* init_partition=NULL);

    ~Partition();

    Partition Optimize();

    Partition Mutate();

    igraph_real_t GetFittness() const;

 private:
    const igraph_t* _graph;
    const int n_nodes;
    const int n_edges;
    igraph_vector_int_t _membership;
    int n_comms; //communities
    igraph_real_t _fittness;
    RandomChooser chooser;

    //https://en.cppreference.com/w/cpp/numeric/random/discrete_distribution/discrete_distribution

    void InitializePartition(double sample_fraction);

    void NewmanSplit(igraph_vector_int_t* indices, igraph_integer_t comm_id_to_split);

    void RandomSplit(igraph_vector_int_t* indices);

    void RandomMerge();
};
