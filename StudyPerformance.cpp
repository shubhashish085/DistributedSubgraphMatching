#include "matchingcommand.h"
#include "BuildTable.h"
#include "graph.h"
#include "backtracking.h"
#include "FilterVertices.h"
#include "GeneratingFilterPlan.h"
#include "Enumeration.h"
#include "ParallelEnumeration.h"
#include "LoadBalancer.h"
#include "wtime.h"
#include <chrono>
#include <limits>
#include <fstream>
#include <mpi.h>


void analyseWorkEstimation(Graph* query_graph, Graph* data_graph){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];
    ui process_count = 2;

    std::cout << "Started Filtering " << std::endl;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];

    std::cout << "Start Vertex : " << start_vertex << std::endl;

    size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    size_t total_workload = 0;

    for(ui i = 0; i < candidates_count[0]; i++){
        total_workload += est_work_array[i];
    }

    std::cout << "Total Estimated Workload : " << total_workload << std::endl;

}


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


void compareBetweenEstimationAndRealCount(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    int world_rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* candidate_limit = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];
    ui process_count = 2;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];    

    size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);
    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    ParallelEnumeration::compareBetweenEstimationAndRealCount(data_graph, query_graph, candidates,
                                                                              candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count);
    end_time = wtime();

    std::cout << "Process : " << world_rank << " - Time " << end_time - start_time << std::endl;
}


void analyseParallelizationWithEvenWorkloadEstimation(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* candidate_limit = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];
    ui process_count = 2;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];    

    size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    size_t* embedding_cnt_array = ParallelEnumeration::exploreWithEvenWorkloadEstimation(data_graph, query_graph, candidates,
                                                                              candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count);
    for(ui idx = 0; idx < process_count; idx++){
        embedding_count += embedding_cnt_array[idx];
    }

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;
}


void analyseParallelizationWithPullBasedLoadBalancing(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

    ui* matching_order = NULL;
    TreeNode* query_tree = NULL;
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* candidate_limit = NULL;
    size_t call_count = 0;
    size_t output_limit = std::numeric_limits<size_t>::max();
    size_t  embedding_count = 0;
    ui* vertex_participating_in_embedding = new ui[data_graph -> getVerticesCount()];
    ui process_count = 2;

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];    

    size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    ParallelEnumeration::exploreWithPullBasedLoadBalancing(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count);

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;
}


//Final Run
/*int main(int argc, char** argv) {

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

    analyseParallelizationWithEvenWorkloadEstimation(query_graph, data_graph, output_performance_file);

    MPI_Finalize();

    double end_time = wtime();
    std::cout << "The time taken is : " << end_time - start_time << std::endl;

}*/


//Test Run
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

    //compareBetweenEstimationAndRealCount(query_graph, data_graph, output_performance_file);
    analyseParallelizationWithPullBasedLoadBalancing(query_graph, data_graph, output_performance_file);
    MPI_Finalize();

    double end_time = wtime();
    std::cout << "The time taken is : " << end_time - start_time << std::endl;

}

