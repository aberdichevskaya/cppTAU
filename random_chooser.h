#pragma once

#include <algorithm>
#include <vector>
#include <random>


class RandomChooser {
 public:
    RandomChooser()
    : gen(std::random_device()) {} 

    template <typename T>
    std::vector<T> RandomChoice(const std::vector<T> &elements, size_t size, 
                            const std::vector<double>& probs = {}) const;

    std::vector<int> RandomChoice(int max_value, size_t size, 
                            const std::vector<double>& probs = {}) const;

 private:
    mutable std::mt19937 gen;
};


template <typename T>
std::vector<T> RandomChooser::RandomChoice(const std::vector<T> &objects, int size, 
                            const std::vector<double>& probs) const {
    if (!probs.empty()) {
        std::vector<T> selected_objects;
        std::vector<size_t> indices(objects.size()); 
        std::vector<T> objects_copy = objects;
        std::iota(indices.begin(), indices.end(), 0);

        for (size_t i = 0; i < size; ++i) {
            std::discrete_distribution<> dist(probs.begin(), probs.end());
            size_t idx = dist(gen); 
            selected_objects.push_back(objects_copy[idx]); 
            
            objects_copy.erase(objects_copy.begin() + idx);
            probs.erase(probs.begin() + idx);
        }
        return selected_objects;
    } else {
        std::vector<T> shuffled = objects;
        std::shuffle(shuffled.begin(), shuffled.end(), gen);
        return std::vector<T>(shuffled.begin(), shuffled.begin() + size);
    }
}

std::vector<int> RandomChooser::RandomChoice(int max_value, int size, 
                            const std::vector<double>& probs) const {
    std::vector<int> range(max_value);
    std::iota(range.begin(), range.end(), 0);
    return RandomChoice(range, size, probs);
}
