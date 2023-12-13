// igraph is a Graph type. igraph/include/igraph_datatype.h.
// igraph_vcount, igraph_ecount
// all includes by #include igraph/include/igraph.h (i think)

// std::discrete_distribution distribution(cbegin(probabilities), cend(probabilities));

int main() {
    std::srand(std::time(0));

    //чтение графов есть в igraph/src/io. Например, igraph_read_graph_edgelist. как быть с тем, что могут быть разные форматы файлов - пока не знаю
    // но в коде Гала вроде как может читаться только список смежности (из любого файла)
    // он использует функцию из nx, а потом уже преобразует это в igraph

    return 0;
}