#pragma once

#include <algorithm>
#include <vector>
#include <random>
#include <stdexcept>

class RandomChooser {
 public:
    RandomChooser()
        : gen(std::random_device{}()) {}  

    template <typename T>
    std::vector<T> RandomChoice(const std::vector<T>& objects, size_t size, 
                                const std::vector<double>& probs = {}) const;

    std::vector<int> RandomChoice(int max_value, size_t size, 
                                  const std::vector<double>& probs = {}) const;

 private:
    mutable std::mt19937 gen;
};

template <typename T>
std::vector<T> RandomChooser::RandomChoice(const std::vector<T>& objects, size_t size, 
                                           const std::vector<double>& probs) const {
    if (objects.size() < size) {
        throw std::length_error("Size requested is larger than the number of objects available.");
    }

    if (size == 0) {
        return {}; 
    }

    if (!probs.empty() && probs.size() != objects.size()) {
        throw std::invalid_argument("Size of probabilities does not match the number of objects.");
    }

    if (!probs.empty()) {
        std::vector<T> selected_objects;
        std::vector<size_t> indices(objects.size()); 
        std::iota(indices.begin(), indices.end(), 0);

        std::vector<T> objects_copy = objects;
        std::vector<double> probs_copy = probs;

        for (size_t i = 0; i < size; ++i) {
            std::discrete_distribution<> dist(probs_copy.begin(), probs_copy.end());
            size_t idx = dist(gen); 
            selected_objects.push_back(objects_copy[idx]); 
            
            objects_copy.erase(objects_copy.begin() + idx);
            probs_copy.erase(probs_copy.begin() + idx);
        }
        return selected_objects;
    } else {
        std::vector<T> shuffled = objects;
        std::shuffle(shuffled.begin(), shuffled.end(), gen);
        return std::vector<T>(shuffled.begin(), shuffled.begin() + size);
    }
}
