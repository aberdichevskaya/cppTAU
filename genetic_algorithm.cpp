#include "genetic_algorithm.h"

//#include <execution>
#include <algorithm>
#include <thread>
#include <iostream>
#include <unordered_map>

// нужно взвесить, насколько мне действительно нужен отдельный RandomChooser

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
                   int32_t n_workers,
                   uint32_t selection_power, 
                   double p_elite, 
                   double p_immigrants, 
                   size_t stopping_criterion_generations,
                   double stopping_criterion_jaccard, 
                   double elite_similarity_threshold)
    : population_size(std::max(population_size, 10)) //почему 10?
    , max_generations(max_generations)
    , n_elite(p_elite*population_size)
    , n_immigrants(n_immigrants*population_size)
    , n_offspring(population_size-n_elite-n_immigrants)
    , stopping_criterion_generations(stopping_criterion_generations)
    , stopping_criterion_jaccard(stopping_criterion_jaccard)
    , elite_similarity_threshold(elite_similarity_threshold) {
        if (n_workers == -1) {
            this->n_workers = std::min(std::thread::hardware_concurrency(), population_size);
        } else {
            this->n_workers = std::min({std::thread::hardware_concurrency(), 
                                        population_size, 
                                        n_workers}, [](uint32_t a, uint32_t b) {
                                            return a < b;
                                        }));
        }
        this->_graph = graph;
        probs = CalculateProbabilities(selection_power, population_size);
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
    std::uniform_real_distribution<double> dis(0.2, 0.9); //вынести это в поля? 
    std::vector<Partition> population(population_size);
    //это должно быть параллельным TODO
    for (size_t i = 0; i < population_size; ++i) {
        population[i] = Partition(_graph, dis(gen));
    }
    return population;
}

void overlap(const igraph_vector_int_t* membership1,  // не сделать ли эту функцию статическим методом, причем наверное даже не geneticalgorithm, а partition? TODO 
                     const igraph_vector_int_t* membership2,
                     igraph_vector_int_t* consensus) {
    // протестировала, вроде работает
    igraph_vector_int_update(consensus, membership1); //точно ли надо делать глубокое копирование. TODO проверить, нужен ли потом ещё membership старого разбиения
    std::unordered_map<igraph_integer_t, igraph_integer_t> consensus_values;
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
    Partition offspring(_graph, 0.5, &partitions_overlap);
    return offspring;
}

std::vector<Partition> GeneticAlgorithm::PopulationCrossover() const {
    //почему некоторые потомки остаются as_is? TO ASK
    std::vector<std::pair<size_t, size_t>> indices_to_cross;
    std::vector<Partition> as_is_offsprings;
    indices_to_cross.reserve(n_offspring);
    as_is_offsprings.reserve(n_offspring);

    for (uint32_t i = 0; i < n_offspring; ++i) {
        auto indices = chooser.RandomChoice(population_size, 2, probs);
        if (FlipCoin()) {
            indices_to_cross.push_back({indices[0], indices[1]});
        } else {
            as_is_offsprings.push_back(_population[indices[0]]); //не нравится мне то, что здесь indices[1] просто игнорируются
            // +может быть стоит мувать, чтобы там membership'ы не копировались
            // так или иначе, лучше у партишенов конструкторы копирования/перемещения написать TODO
        }
        std::vector<Partition> crossed_offsprings(indices_to_cross.size());
        //TODO распараллелить
        for (size_t i = 0; i < indices_to_cross.size(); ++i) {
            crossed_offsprings[i] = SingleCrossover(indices_to_cross[i].first, 
                                                    indices_to_cross[i].second);
        }
        crossed_offsprings.insert(crossed_offsprings.back(), as_is_offsprings.begin(),
                                    as_is_offsprings.end()); //TODO std::make_move_iterator?
        // мне не то чтобы нравится такое переиспользование crossed_offspring, но без move пока так
        return crossed_offsprings;
    }
}


// Жаккар проще для понимания, О(1) памяти. Ассимптотика О(n^2)
// ARI сложнее, O(k^2), ассимптотика O(n^2), но оптимизируется до O(nk) или O(n). Точнее и чувствительнее
// TODO сравнить их на практике и выбрать оптимальный
// либо добавить пользователю возможность выбирать, но это чёт сложно
// TODO оба индекса хорошо параллелятся. Хотя в этом не особо есть смысл, это будет вложенная параллелизация
inline double comb2(int n) {
    return static_cast<double>(n) * (n - 1) / 2.0;
}

// Функция для вычисления ARI 
double ComputePartitionSimilarityARI(const igraph_vector_int_t* membership1, 
                                const igraph_vector_int_t* membership2) { //static метод Partition или GeneticAlgorithm? TODO
    size_t n = igraph_vector_int_size(membership1);
    std::unordered_map<igraph_integer_t, int32_t> count_membership1, count_membership2;
    std::unordered_map<std::pair<igraph_integer_t, igraph_integer_t>, int32_t> count_joint;

    for (size_t i = 0; i < n; ++i) {
        igraph_integer_t membership1_i = membership1->stor_begin[i];
        igraph_integer_t membership2_i = membership2->stor_begin[i];
        count_membership1[membership1_i]++;
        count_membership2[membership2_i]++;
        count_joint[{membership1_i, membership2_i}]++;
    }

    double sum_comb_memb1 = 0, sum_comb_memb2 = 0, sum_comb_joint = 0;
    for (const auto& el : count_membership1) {
        sum_comb_memb1 += comb2(el.second);
    }
    for (const auto& el : count_membership2) {
        sum_comb_memb2 += comb2(el.second);
    }
    for (const auto& el : count_joint) { //вот тут может вылезти O(n^2)
        sum_comb_joint += comb2(el.second);
    }

    double index = (sum_comb_joint - (sum_comb_memb1 * sum_comb_memb2) / comb2(n)) / ((sum_comb_memb1 + sum_comb_memb2) / 2 - (sum_comb_memb1 * sum_comb_memb2) / comb2(n));
    return index;
}

double ComputePartitionSimilarityJaccard(const igraph_vector_int_t* membership1, 
                                const igraph_vector_int_t* membership2) {
    int a = 0, c = 0, d = 0;
    size_t n = igraph_vector_int_size(membership1);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            bool in_same_cluster_1 = membership1->stor_begin[i] == membership1->stor_begin[j];
            bool in_same_cluster_2 = membership2->stor_begin[i] == membership2->stor_begin[j];
            if (in_same_cluster_1 && in_same_cluster_2) {
                a++;
            } else if (in_same_cluster_1 || in_same_cluster_2) {
                c++;
            } else {
                d++;
            }
        }
    }
    double jac = static_cast<double>(a) / (a + c + d);
    return jac;
}


std::vector<size_t> GeneticAlgorithm::EliteSelection() const {
// пока вообще без оптимизаций
// TODO возможные оптимизации: параллельность; проверка, что не происходит повторное пересчитывание similarity
// для одной и той же пары разбиений; ограничение на количество итераций (и опция докинуть просто случайные разбиения в элиту)
    std::vector<size_t> elite_indices = {0}; // starting with an assumption that [0] is already an elite (because population is sorted)
    size_t candidate_idx = 1; // element index to start checking
    while (elite_indices.size() < n_elite && candidate_idx < _population.size()) {
        bool is_elite = true;
        for (auto elite_idx : elite_indices) {
            double similarity =  ComputePartitionSimilarityJaccard( // ComputePartitionSimilarityARI(
                                    _population[elite_idx].GetMembership(),
                                    _population[candidate_idx].GetMembership()); 
            if (similarity > elite_similarity_threshold) {
                is_elite = false;
                break;
            }
        }
        if (is_elite) {
            elite_indices.push_back(candidate_idx);
        }
        candidate_idx++;
    }
    return elite_idx;
}

std::pair<Partition, std::vector<double>> GeneticAlgorithm::Run() {
    std::vector<double> best_modularity_per_generation; //неплохо бы сколько-то зарезервировать, но сколько?
    int32_t cnt_convergence = 0;
    std::vector<igraph_integer_t> last_best_partition;
    last_best_partition.reserve(igraph_vcount(_graph));
    
    _population = CreatePopulation(population_size);
    for (int32_t generation_i = 1; generation_i <= max_generations; ++generation_i) {
        // TODO добавить замер времени
        for (auto& indiv : _population) { // TODO распараллелить
            indiv.Optimize();
        }

        std::stable_sort(/*std::execution::parallel_policy,*/ _population.begin(), _population.end(),
                            [](const Partition &a, const Partition &b) {
                                return a.GetFittness() > b.GetFittness();
                            });
        // не то чтобы параллельная сортировка очень нужна, популяция маленькая, так что надо проверять TODO

        Partition best_indiv = _population[0];
        double best_score = best_indiv.GetFittness();

        if (last_best_memb.size()) {
            double sim_to_last_best = ComputePartitionSimilarityJaccard( // ComputePartitionSimilarityARI(
                                    best_indiv.GetMembership(), last_best_memb);
            if (sim_to_last_best > stopping_criterion_jaccard) {
                cnt_convergence++;
            } else {
                cnt_convergence = 0;
            }
        }
        last_best_memb = best_indiv.GetMembership();
        if (cnt_convergence == stopping_criterion_generations || generation_i == max_generations) {
            break;
        }

        std::vector<size_t> elite_indices = EliteSelection();
        std::vector<Partition> elite(n_elite);
        for (size_t idx : elite_indices) {
            elite[i] = _population[idx];
        }

        // TODO parallel
        std::vector<Partition> offsprings = PopulationCrossover();
        for (auto& offspring : offsprings) {
            offspring.Mutate();
        }

        std::vector<Partition> immigrants = CreatePopulation(n_immigrants);

        // опять же, параллельность не то чтобы очень оправдана, но компилятор должен уметь самостоятельно принять это решение
        std::copy(/*std::execution::parallel_policy, */std::make_move_iterator(elite.begin()), 
                                std::make_move_iterator(elite.end()), _population.begin());
        std::copy(/*std::execution::parallel_policy, */std::make_move_iterator(immigrants.begin()), 
            std::make_move_iterator(immigrants.end()), _population.begin() + n_elite);
        std::copy(/*std::execution::parallel_policy, */std::make_move_iterator(offsprings.begin()), 
            std::make_move_iterator(offsprings.end()), _population.begin() + n_elite + n_immigrants);
        
        std::cout << "Generation " << generation_i << ". Best score: " << best_score << "\n";
    }
}
