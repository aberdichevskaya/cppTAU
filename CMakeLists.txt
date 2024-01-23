cmake_minimum_required(VERSION 3.18)

set(PROJECT_NAME cppTAU)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP REQUIRED)

option(USE_SANITIZER "Use address sanitizer" OFF)
option(USE_PROFILER "Use profiler" OFF)

if(USE_SANITIZER MATCHES address AND CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Using Address Sanitizer")
    add_compile_options(-fsanitize=address)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -fsanitize=address")
elseif(USE_SANITIZER MATCHES thread AND CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Using Thread Sanitizer")
    add_compile_options(-fsanitize=thread)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g  -fsanitize=thread")
elseif(USE_PROFILER MATCHES ON AND CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Using Profiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
    message(STATUS "Debug build type: Adding -g flag for debugging")
    add_compile_options(-g)
else()
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
