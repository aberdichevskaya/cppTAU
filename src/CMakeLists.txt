add_subdirectory(genetic_algorithm)
add_subdirectory(partition)
add_subdirectory(random_chooser)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE igraph genetic_algorithm)

target_include_directories(${PROJECT_NAME} PRIVATE
                            ${CMAKE_SOURCE_DIR}/vendor/argparse/
                            ${CMAKE_SOURCE_DIR}/vendor/igraph/include/
                        )
