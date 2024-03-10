set(PARTITION_LIBRARY partition)

message("${CMAKE_SOURCE_DIR}")

add_library(${PARTITION_LIBRARY} STATIC partition.cpp)

target_link_libraries(${PARTITION_LIBRARY} PRIVATE igraph random_chooser)

target_include_directories(${PARTITION_LIBRARY} PRIVATE 
                ${CMAKE_SOURCE_DIR}/vendor/igraph/include/)
