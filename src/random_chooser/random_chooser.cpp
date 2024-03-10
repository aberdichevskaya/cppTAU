#include "random_chooser.h"

std::vector<int> RandomChooser::RandomChoice(int max_value, size_t size, 
                                             const std::vector<double>& probs) const {
    if (max_value <= 0) {
        throw std::invalid_argument("max_value must be positive.");
    }
    if (size_t(max_value) < size) {
        throw std::length_error("Size requested is larger than the max_value.");
    }
    std::vector<int> range(max_value);
    std::iota(range.begin(), range.end(), 0);
    return RandomChoice(range, size, probs);
}
