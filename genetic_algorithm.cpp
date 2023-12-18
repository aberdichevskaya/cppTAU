#include "genetic_algorithm.h"

#include <execution>
#include <algorithms>
#include <thread>
#include <iostream>

GeneticAlgorithm::GeneticAlgorithm(igraph_t* graph, 
                   size_t population_size,
                   size_t max_generations, 
                   size_t n_workers,
                   uint32_t selection_power, 
                   double p_elite, 
                   double p_immigrants, 
                   size_t stopping_criterion_generations,
                   double stopping_criterion_jaccard, 
                   double elite_similarity_threshold)
    : population_size(max(population_size, 10)) //почему 10?
    , max_generations(max_generations)
    , n_workers(min({std::thread::hardware_concurrency(), population_size, n_workers})) 
    , n_elite(p_elite*population_size)
    , n_immigrants(n_immigrants*population_size)
    , stopping_criterion_generations(stopping_criterion_generations)
    , stopping_criterion_jaccard(stopping_criterion_jaccard)
    , elite_similarity_threshold(elite_similarity_threshold) {
        this->_graph = graph;
        CalculateProbabilities(selection_power);
        std::cout << "Main parameter values: pop_size = " << population_size << 
            ", workers = " << n_workers << ", max_generations = " << max_generations << "\n";
}


int64_t binpow(int64_t x, int32_t a) {
    int64_t res = 1;
    for (; a > 0; a >>= 1) {
        if (a & 1) {
            res *= x;
        }
        res *= res;
    }
    return res;
}

void GeneticAlgorithm::CalculateProbabilities(uint32_t selection_power) {
    probs.resize(population_size);
    std::vector<int64_t> values(population_size);
    std::iota(values.rbegin(), values.rend(), 1);
    int64_t denom = 0;
    for (auto& value : values) {
        value = binpow(value, selection_power);
        denom += value;
    }
    for (int i = 0; i < population_size; ++i) {
        probs[i] = static_cast<double>(values[i]) / denom;
    }
}

std::vector<Partition> CreatePopulation(size_t population_size) const {
    
}
