#pragma once

#include "igraph/include/igraph.h"

#include "partition.h"
#include "random_chooser.h"

#include <vector>
#include <random>

class GeneticAlgorithm {
 public:
    GeneticAlgorithm(const igraph_t* graph, 
                   size_t population_size,
                   size_t max_generations,
                   size_t n_workers, 
                   uint32_t selection_power=5, 
                   double p_elite=0.1, 
                   double p_immigrants=0.15, 
                   size_t stopping_criterion_generations=10,
                   double stopping_criterion_jaccard=0.98, 
                   double elite_similarity_threshold=0.9);

    std::pair<Partition, std::vector<double>> Run(); /*альтернативное название - FindPartition*/

 private:
    void CreatePopulation(std::vector<Partition> &population) const;
    Partition SingleCrossover(size_t idx1, size_t idx2) const;
    std::vector<Partition> PopulationCrossover() const;
    std::vector<size_t> EliteSelection() const;

    const igraph_t* _graph;
    const size_t population_size;
    const size_t max_generations; 
    size_t n_workers;
    const size_t n_elite;
    const size_t n_immigrants;
    const size_t n_offspring;
    const size_t stopping_criterion_generations;
    const double stopping_criterion_jaccard;
    const double elite_similarity_threshold;

    std::vector<Partition> _population; 

    RandomChooser chooser;

    std::vector<double> probs;
};
