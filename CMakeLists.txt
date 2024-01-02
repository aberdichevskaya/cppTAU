cmake_minimum_required(VERSION 3.18)

set(PROJECT_NAME cppTAU)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

option(USE_SANITIZER "Use address sanitizer" OFF)

if(USE_SANITIZER MATCHES address)
    message(STATUS "Using Address Sanitizer")
    add_compile_options(-g -fsanitize=address)
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug build type: Adding -g flag for debugging")
    add_compile_options(-g)
elseif(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Release)
    message(STATUS "Setting default build type to Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Set build type to default (Release)." FORCE)
endif()

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

add_subdirectory(igraph)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC igraph/)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC igraph)

# Если igraph зависит от других библиотек, вам также может потребоваться их подключить.
# Например:
# target_link_libraries(${PROJECT_NAME} m gmp glpk)
