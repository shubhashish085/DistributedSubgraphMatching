#include "matchingcommand.h"
#include "BuildTable.h"
#include "graph.h"
#include "backtracking.h"
#include "FilterVertices.h"
#include "GeneratingFilterPlan.h"
#include "Enumeration.h"
#include "ParallelEnumeration.h"
#include "wtime.h"
#include <chrono>
#include <limits>
#include <fstream>
#include <mpi.h>


void analyseParallelization(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];
    ui process_count = 2;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    std::cout << "####### Candidate count  : " ;

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << candidates_count[i] << " " ;
    }

    std::cout << std::endl;

    VertexID start_vertex = matching_order[0];

    ui* cand_degree_offset = new ui[candidates_count[start_vertex] + 1];
    cand_degree_offset[0] = 0;
    ui* candidate_limit = NULL;

    for(ui j = 1; j < candidates_count[start_vertex] + 1; j++){
        cand_degree_offset[j] = cand_degree_offset[j - 1] + data_graph->getVertexDegree(candidates[start_vertex][j - 1]);
    }

    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    size_t* embedding_cnt_array = ParallelEnumeration::exploreWithEvenDegreeDist(data_graph, query_graph, candidates,
                                                                              candidates_count, matching_order, query_tree, output_limit, call_count);
    for(ui idx = 0; idx < process_count; idx++){
        embedding_count += embedding_cnt_array[idx];
    }

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;


}

int main(int argc, char** argv) {

    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string output_performance_file = command.getOutputFilePath();


    std::cout << " Query Graph : " << input_query_graph_file << std::endl;
    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->printGraphMetaData();

    std::cout << " Data Graph : " << input_data_graph_file << std::endl;
    Graph* data_graph = new Graph();
    data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    data_graph->printGraphMetaData();

    double start_time = wtime();

    MPI_Init(NULL, NULL);

    analyseParallelization(query_graph, data_graph, output_performance_file);

    MPI_Finalize();

    double end_time = wtime();
    std::cout << "The time taken is : " << end_time - start_time << std::endl;

}


