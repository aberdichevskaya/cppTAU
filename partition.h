#pragma once

#include "random_chooser.h"
#include "igraph/include/igraph.h"

#include <ctime>

#include <vector>
#include <iostream>
#include <mutex>

class Partition {
 public:
    Partition(const igraph_t* graph, double sample_fraction,
                    const igraph_vector_int_t* init_partition=NULL);
    
    Partition();

    Partition(const Partition& p);

    Partition& operator=(const Partition& p);
    
    ~Partition();

    void Optimize();

    void Mutate();

    igraph_real_t GetFittness() const {
      std::lock_guard<std::mutex> lock(mtx);
      return _fittness;
    }

    const igraph_vector_int_t* GetMembership() const { // лочить внутри этого метода мьютекс смысла нет. надо какую-то оболочку использовать для вектора или хз что
      //возможно, стоит вообще не выдавать доступ к _membership (убрать этот метод), и сделать внутренние методы, которые будут таскать _membership + какие-то методы, которые сами обращаются к индексу в _membership
      return &_membership;
    }

 private:
    const igraph_t* _graph;
    size_t n_nodes;
    size_t n_edges;
    igraph_vector_int_t _membership;
    igraph_integer_t n_comms; //communities
    igraph_real_t _fittness;
    RandomChooser chooser;
    mutable std::mutex mtx;

    void InitializePartition(double sample_fraction);

    void NewmanSplit(igraph_vector_int_t* indices, igraph_integer_t comm_id_to_split);

    void RandomSplit(igraph_vector_int_t* indices);

    void RandomMerge(size_t edges_subsample_size);

    void CheckMembershipSize();
};
