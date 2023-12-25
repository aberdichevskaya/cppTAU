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

    //чтение графов есть в igraph/src/io. Например, igraph_read_graph_edgelist. как быть с тем, что могут быть разные форматы файлов - пока не знаю
    // но в коде Гала вроде как может читаться только список смежности (из любого файла)
    // он использует функцию из nx, а потом уже преобразует это в igraph

    argparse::ArgumentParser program("cppTAU");
    program.add_argument("--graph")
        .required()
        .help("path to graph file; supports adjacency list format.");
    program.add_argument("--size")
        .default_value(60)
        .help("size of population; default is 60.");
    program.add_argument("--workers")
        .default_value(-1) //это тупа какое-то слишком большое значение, мб можно поэлегантнее придумать
        .help("number of workers; default is number of available CPUs.");
    program.add_argument("--max_generations")
        .default_value(500)
        .help("maximum number of generations to run; default is 500.");

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        std::exit(1);
    }

    auto file_path = program.get<char*>("--graph");
    FILE *file = fopen(file_path, "r");
    igraph_t graph;
    igraph_read_graph_edgelist(&graph, file, 0, false);
    fclose(file);

    GeneticAlgorithm algorithm(&graph, 
                   program.get<size_t>("--size"),
                   program.get<size_t>("--max_generations"),
                   program.get<int32_t>("--workers"));
    auto res = algorithm.Run();



    return 0;
}