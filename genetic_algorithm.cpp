#include "genetic_algorithm.h"

#include <execution>
#include <algorithms>
#include <thread>
#include <iostream>
#include <unordered_map>

std::vector<double> CalculateProbabilities(uint32_t selection_power, size_t population_size) {
    std::vector<double> probs(population_size);
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
    return probs;
}

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
    , n_offspring(population_size-n_elite-n_immigrants)
    , stopping_criterion_generations(stopping_criterion_generations)
    , stopping_criterion_jaccard(stopping_criterion_jaccard)
    , elite_similarity_threshold(elite_similarity_threshold) {
        this->_graph = graph;
        auto probs = CalculateProbabilities(selection_power, population_size);
        dis = std::discrete_distribution<size_t>(probs.begin(), probs.end());
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

std::vector<Partition> GeneticAlgorithm::CreatePopulation(size_t population_size) const {
    std::mt19937 gen(std::random_device());
    std::uniform_real_distribution<double> dis(0.2, 0.9); //вынести это в поля? или сделать метод статическим?
    std::vector<Partition> population(population_size);
    //это должно быть параллельным TODO
    for (size_t i = 0; i < population_size; ++i) {
        population[i] = Partition(_graph, dis(gen));
    }
    return population;
}

void overlap(const igraph_vector_int_t* membership1, 
                     const igraph_vector_int_t* membership2,
                     igraph_vector_int_t* consensus) {
    // протестировала, вроде работает
    igraph_vector_int_update(consensus, membership1); //точно ли надо делать глубокое копирование. TODO проверить, нужен ли потом ещё membership старого разбиения
    unordered_map<igraph_integer_t, igraph_integer_t> consensus_values;
    igraph_integer_t comm = 0;
    // TODO распараллелить
    for (size_t i = 0; i < igraph_vector_int_size(membership1); ++i) {
        if (membership1->stor_begin[i] == membership2->stor_begin[i]) {
            if (!consensus_values.count(membership1->stor_begin[i])) {
                consensus_values[membership1->stor_begin[i]] = comm++;  
            } 
            consensus->stor_begin[i] = consensus_values[membership1->stor_begin[i]];
        } else {
            consensus->stor_begin[i] = comm++;
        }
    }
}

Partition GeneticAlgorithm::SingleCrossover(size_t idx1, size_t idx2) const {
    auto membership1 = _population[idx1].GetMembership();
    auto membership2 = _population[idx2].GetMembership();
    igraph_vector_int_t partitions_overlap;
    igraph_vector_int_init(&partitions_overlap, igraph_vector_int_size(membership1));
    overlap(membership1, membership2, &partitions_overlap);
    Partition offspring(_graph, 0.5, partitions_overlap);
    return offspring;
}

std::vector<Partition> GeneticAlgorithm::PopulationCrossover() const {
    //почему некоторые потомки остаются as_is? TO ASK
    std::vector<std::pair<size_t, size_t>> indices_to_cross;
    std::vector<Partition> as_is_offspring;
    indices_to_cross.reserve(n_offspring);
    as_is_offsprings.reserve(n_offspring);

    for (uint32_t i = 0; i < n_offspring; ++i) {
        auto indices = chooser.RandomChoice(population_size, 2, dis);
        if (FlipCoin()) {
            indices_to_cross.push_back({indices[0], indices[1]});
        } else {
            as_is_offsprings.push_back(population[indices[0]]); //не нравится мне то, что здесь indices[1] просто игнорируются
            // +может быть стоит мувать, чтобы там membership'ы не копировались
            // так или иначе, лучше у партишенов конструкторы копирования/перемещения написать TODO
        }
        std::vector<Partition> crossed_offsprings(indices_to_cross.size());
        //TODO распараллелить
        for (size_t i = 0; i < indices_to_cross.size(); ++i) {
            crossed_offspring[i] = SingleCrossover(indices_to_cross[i].first, 
                                                    indices_to_cross[i].second);
        }
        crossed_offspring.insert(crossed_offspring.back(), as_is_offspring.begin(),
                                    as_is_offspring.end()); //TODO move_iterator?
        // мне не то чтобы нравится такое переиспользование crossed_offspring, но без move пока так
        return crossed_offspring;
    }
}

// мне не кажется адекватной идея использовать расстояние жаккара. chatgpt считает, что лучше брать
// метрики ARI и NMI. нужно попытаться понять, почему Гал решил использовать Жаккара, и точно
// ли эти 2 лучше
