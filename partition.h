#pragma once

#include "random_chooser.h"
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
    Partition(const igraph_t* graph, double sample_fraction,
                    const igraph_vector_int_t* init_partition=NULL);
    
    Partition();

    Partition& operator=(const Partition& p);
    
    ~Partition();

    void Optimize();

    void Mutate();

    igraph_real_t GetFittness() const {
      return _fittness;
    }

    const igraph_vector_int_t* GetMembership() const { // точно const?
      return &_membership;
    }

 private:
    const igraph_t* _graph;
    size_t n_nodes;
    size_t n_edges;
    igraph_vector_int_t _membership;
    bool membership_is_inialized;
    igraph_integer_t n_comms; //communities
    igraph_real_t _fittness;
    RandomChooser chooser;

    void InitializePartition(double sample_fraction);

    void NewmanSplit(igraph_vector_int_t* indices, igraph_integer_t comm_id_to_split);

    void RandomSplit(igraph_vector_int_t* indices);

    void RandomMerge(size_t edges_subsample_size);
};
