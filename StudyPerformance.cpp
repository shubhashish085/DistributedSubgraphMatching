#include "matchingcommand.h"
#include "BuildTable.h"
#include "graph.h"
#include "backtracking.h"
#include "FilterVertices.h"
#include "Automorphism.h"
#include "GeneratingFilterPlan.h"
#include "Enumeration.h"
#include "ParallelEnumeration.h"
#include "LoadBalancer.h"
#include "wtime.h"
#include "GapMeasure.h"
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


void analyseParallelizationWithPushBasedLoadBalancing(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

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
    ParallelEnumeration::exploreWithPushBasedLoadBalancing(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count);

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;
}


void analysePushBasedLoadBalancingWithNoWaiting(Graph* query_graph, Graph* data_graph, const std::string& output_file_path){

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
    ParallelEnumeration::explorePushBasedLoadBalancingWithNoWaiting(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count);

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;
}

void print_schedule_restriction_map(std::map<ui, std::vector<std::pair<ui, ui>>>& schedule_restriction_map, ui size){

    for(int i = 0; i < size; i++){
        std::map<ui, std::vector<std::pair<ui, ui>>>::iterator it = schedule_restriction_map.find(i);

        if(it != schedule_restriction_map.end()){
            std::vector<std::pair<ui, ui>>::iterator vtr_itr = (it->second).begin();
            std::cout << "Schedule at index : " << i << " : ";

            while (vtr_itr != (it->second).end()){
                std::cout << "(" << vtr_itr -> first << ", " <<  vtr_itr -> second << ") ";
                ++vtr_itr;
            }

            std::cout << std::endl;
        }
    }
}

void analyseAutomorphismBreak(Graph* query_graph, Graph* data_graph){

    std::vector< std::pair<ui, ui> > ordered_pairs;
    std::map<ui, std::vector<std::pair<ui, ui>>> schedule_restriction_map; 
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

    ui* adj_mat = Automorphism::convert_to_adj_mat(query_graph-> getVerticesCount(), query_graph->getOffsets(), query_graph ->getNeighbors());

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];

    Automorphism::aggressive_optimize(ordered_pairs, adj_mat, query_graph->getVerticesCount());
    Automorphism::restriction_integration_with_scheduling(matching_order, query_graph->getVerticesCount(), ordered_pairs, schedule_restriction_map); 

    print_schedule_restriction_map(schedule_restriction_map, query_graph->getVerticesCount());  

    size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    ParallelEnumeration::exploreGraphWithAutomorphismBreak(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree, est_work_array, output_limit, call_count, schedule_restriction_map);

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;
}

void analyseHybridParallelization(Graph* query_graph, std::string data_graph_file){


    int world_rank, allocated_file_number;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(world_rank == 0){
        allocated_file_number = 0;
    }else{
        allocated_file_number = world_rank - 1;
    }


    std::string input_data_file = data_graph_file + "_" + std::to_string(allocated_file_number) + ".graph";
    Graph* data_graph = new Graph();
    data_graph->loadGraphFromFileWithReindexing(input_data_file);
    data_graph->printGraphMetaData();
    
    
    std::vector< std::pair<ui, ui> > ordered_pairs;
    std::map<ui, std::vector<std::pair<ui, ui>>> schedule_restriction_map; 
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

    ui* adj_mat = Automorphism::convert_to_adj_mat(query_graph-> getVerticesCount(), query_graph->getOffsets(), query_graph ->getNeighbors());

    FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    VertexID start_vertex = matching_order[0];

    //Automorphism::aggressive_optimize(ordered_pairs, adj_mat, query_graph->getVerticesCount());
    //Automorphism::restriction_integration_with_scheduling(matching_order, query_graph->getVerticesCount(), ordered_pairs, schedule_restriction_map); 

    //print_schedule_restriction_map(schedule_restriction_map, query_graph->getVerticesCount());  

    //size_t* est_work_array = LoadBalancer::workloadEstimator(data_graph, query_graph, candidates, candidates_count, matching_order, query_tree);

    //Parallel Strategy
    double start_time, end_time;

    
    embedding_count = 0;
    call_count = 0;

    start_time = wtime();
    ParallelEnumeration::exploreGraphInHybridFashion(data_graph, query_graph, candidates, candidates_count, matching_order,
                                                        query_tree, output_limit, call_count, schedule_restriction_map);

    end_time = wtime();

    std::cout << "Time " << end_time - start_time << std::endl;

}


//Main
// int main(int argc, char** argv) {

//     MatchingCommand command(argc, argv);
//     std::string input_query_graph_file = command.getQueryGraphFilePath();
//     std::string input_data_graph_file = command.getDataGraphFilePath();
//     std::string output_performance_file = command.getOutputFilePath();


//     std::cout << " Query Graph : " << input_query_graph_file << std::endl;
//     Graph* query_graph = new Graph();
//     query_graph->loadGraphFromFile(input_query_graph_file);
//     query_graph->printGraphMetaData();

//     std::cout << " Data Graph : " << input_data_graph_file << std::endl;
//     Graph* data_graph = new Graph();
//     data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
//     data_graph->printGraphMetaData();

//     double start_time = wtime();

//     MPI_Init(NULL, NULL);

//     analyseAutomorphismBreak(query_graph, data_graph);
//     MPI_Finalize();

//     double end_time = wtime();
//     std::cout << "The time taken is : " << end_time - start_time << std::endl;

// }

//For Partitioned Graph
/*int main(int argc, char** argv) {

    std::string prefix_file_name = "orkut";
    std::string input_data_graph_directory = "/home/kars1/Research_Projects/metis/";
    std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/com-orkut.ungraph.txt";
    

    int numberOfParts[3] = {2, 4, 8};

    
    MPI_Init(NULL, NULL);

    Graph* main_data_graph = new Graph();
    main_data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    main_data_graph->printGraphMetaData();


    for(int i = 0; i < 3; i++){

        long long gap_distance = 0;
        double epsilon = 0;

        for(int j = 0; j < numberOfParts[i]; j++){
            std::string filename = input_data_graph_directory + std::to_string(numberOfParts[i]) + "_" + prefix_file_name + "_" +  std::to_string(j) + ".graph";
            std::cout << " Data Graph : " << filename << std::endl;
            Graph* data_graph = new Graph();
            data_graph->loadGraphFromFileWithReindexing(filename);
            data_graph->printGraphMetaData();

            gap_distance += GapMeasure::measure_the_epsilon_gap(data_graph);
        }

        std::string filename = input_data_graph_directory + std::to_string(numberOfParts[i]) + "_" + prefix_file_name + "_partition" + ".graph";
        std::cout << " Data Graph : " << filename << std::endl;
        Graph* data_graph = new Graph();
        gap_distance += data_graph->measureGapForPartitionedEdges(filename);
        epsilon = (double)(1.0 * gap_distance) / (main_data_graph -> getEdgesCount());  

        std::cout << "Graph : " << numberOfParts[i] << prefix_file_name << std::endl;
        std::cout << "Epsilon : " <<  epsilon << std::endl;           

    }

    MPI_Finalize();    

}*/


//Vertex Reordering for the entire graph
int main(int argc, char** argv) {

    std::string input_data_graph_file = "/home/kars1/Parallel_computation/dataset/com-orkut.ungraph.txt";

    
    MPI_Init(NULL, NULL);

    Graph* main_data_graph = new Graph();
    main_data_graph->loadGraphFromFileWithoutStringConversion(input_data_graph_file);
    main_data_graph->printGraphMetaData();
    long long gap_distance = GapMeasure::measure_the_beta_gap(main_data_graph);
    double epsilon = (double)(1.0 * gap_distance) / (main_data_graph -> getVerticesCount());

    std::cout << "Graph : " << input_data_graph_file << std::endl;
    std::cout << "Beta : " <<  epsilon << std::endl;           


    MPI_Finalize();    

}



// Partition Run
/*int main(int argc, char** argv) {

    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_file = command.getDataGraphFilePath();
    std::string output_performance_file = command.getOutputFilePath();


    std::cout << " Query Graph : " << input_query_graph_file << std::endl;
    Graph* query_graph = new Graph();
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->printGraphMetaData();


    double start_time = wtime();

    MPI_Init(NULL, NULL);

    analyseHybridParallelization(query_graph, input_data_graph_file);
    MPI_Finalize();

    double end_time = wtime();
    std::cout << "The time taken is : " << end_time - start_time << std::endl;

}*/

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

    //compareBetweenEstimationAndRealCount(query_graph, data_graph, output_performance_file);
    //analyseParallelizationWithPushBasedLoadBalancing(query_graph, data_graph, output_performance_file);
    analysePushBasedLoadBalancingWithNoWaiting(query_graph, data_graph, output_performance_file);
    MPI_Finalize();

    double end_time = wtime();
    std::cout << "The time taken is : " << end_time - start_time << std::endl;

}*/

