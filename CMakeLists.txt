cmake_minimum_required(VERSION 2.8)
project(DistributedSubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

find_package(MPI REQUIRED)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3 -g -Wall -march=native -pthread")

set(DISTRIBUTED_SUBGRAPH_MATCHING_FILE_LIST
        backtracking.cpp
        backtracking.h
        Automorphism.h
        Automorphism.cpp
        BuildTable.h
        BuildTable.cpp
        FilterVertices.h
        FilterVertices.cpp
        commandparser.h
        commandparser.cpp
        matchingcommand.h
        matchingcommand.cpp
        GeneratingFilterPlan.h
        GeneratingFilterPlan.cpp
        Enumeration.h
        Enumeration.cpp
        ParallelEnumeration.h
        ParallelEnumeration.cpp
        LoadBalancer.h
        LoadBalancer.cpp
        graph.cpp
        graph.h
        utilities.cpp
        utilities.h
        types.h
        util.h
        wtime.h
        StudyPerformance.cpp)


add_executable(SubgraphMatching.out ${DISTRIBUTED_SUBGRAPH_MATCHING_FILE_LIST})

target_include_directories(SubgraphMatching.out PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(SubgraphMatching.out PUBLIC ${MPI_CXX_LIBRARIES})
