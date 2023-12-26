cmake_minimum_required(VERSION 3.18)
project(cppTAU)

set(CMAKE_CXX_STANDARD 17)


set(SOURCES
    main.cpp
    partition.cpp
    genetic_algorithm.cpp
    random_chooser.cpp
)


set(HEADERS
    random_chooser.h
    partition.h
    genetic_algorithm.h
    argparse.hpp
)


include_directories(igraph/) 


add_subdirectory(igraph)

add_executable(cppTAU ${SOURCES} ${HEADERS})

target_link_libraries(cppTAU igraph)

# Если igraph зависит от других библиотек, вам также может потребоваться их подключить.
# Например:
# target_link_libraries(cppTAU m gmp glpk)
