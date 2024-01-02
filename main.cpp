#include "genetic_algorithm.h"
#include "igraph/include/igraph.h"
#include "argparse.hpp"
#include "partition.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>

#include <stdio.h>

int main(int argc, char* argv[]) {
    std::srand(std::time(0));

    argparse::ArgumentParser program("cppTAU");
    program.add_argument("--graph")
        .required()
        .help("path to graph file; supports adjacency list format.");
    program.add_argument("--size")
        .default_value(60)
        .scan<'i', int>()
        .help("size of population; default is 60.");
    program.add_argument("--workers")
        .default_value(0) 
        .scan<'i', int>()
        .help("number of workers; default is number of available CPUs.");
    program.add_argument("--max_generations")
        .default_value(500)
        .scan<'i', int>()
        .help("maximum number of generations to run; default is 500.");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto file_path = program.get<std::string>("--graph");
    FILE *file = fopen(file_path.c_str(), "r");
    igraph_t graph;
    igraph_read_graph_edgelist(&graph, file, 0, false); //граф где-то утекает 
    fclose(file);

    GeneticAlgorithm algorithm(&graph, 
                   program.get<int>("--size"),
                   program.get<int>("--max_generations"),
                   program.get<int>("--workers"));
    auto res = algorithm.Run();
    igraph_destroy(&graph);

    return 0;
}