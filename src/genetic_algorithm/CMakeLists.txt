set(GENETIC_ALGORITHM_LIBRARY genetic_algorithm)

add_library(${GENETIC_ALGORITHM_LIBRARY} STATIC genetic_algorithm.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${GENETIC_ALGORITHM_LIBRARY} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${GENETIC_ALGORITHM_LIBRARY} PRIVATE 
                        igraph partition random_chooser)

target_include_directories(${GENETIC_ALGORITHM_LIBRARY} PRIVATE 
                        ${CMAKE_SOURCE_DIR}/vendor/igraph/include/)
                        