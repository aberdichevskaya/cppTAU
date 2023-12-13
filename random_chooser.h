// не протестировано
// choice может быть const? скорее нет, потому что генератор передаётся в shuffle по &&

#pragma once

#include <algoritm>
#include <vector>
#include <random>


class RandomChooser {
 public:
    RandomChooser()
    : gen(std::time(0)) {} 

    template <typename T, typename Distribution = std::uniform_int_distribution<std::size_t>>
    std::vector<T> RandomChoice(const std::vector<T> &elements, int size, 
                            Distribution distribution=Distribution{});

    template <typename Distribution = std::uniform_int_distribution<std::size_t>>
    std::vector<int> RandomChoice(int max_value, int size, 
                            Distribution distribution=Distribution{});

 private:
    mutable std::mt19937 gen;
};


template <typename T, typename Distribution = std::uniform_int_distribution<std::size_t>>
std::vector<T> RandomChooser::RandomChoice(const std::vector<T> &elements, int size, 
                            Distribution distribution) {
    std::vector<T> shuffled = elements;
    std::shuffle(shuffled.begin(), shuffled.end(), 
        [&distribution, &gen] {
            return distribution(gen);
        });
    return {shuffled.begin(), shuffled.begin() + size};
}


template <typename Distribution = std::uniform_int_distribution<std::size_t>>
    std::vector<int> RandomChoice(int max_value, int size, 
                            Distribution distribution) {
    std::vector<int> range(size);
    std::iota(range.begin(), range().end(), 0);
    return RandomChoice(range, size, distribution);
}
