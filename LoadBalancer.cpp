
#include <fstream>
#include "types.h"
#include "LoadBalancer.h"



size_t LoadBalancer::calculateWorkLoad(const Graph *data_graph, const Graph *query_graph, ui max_valid_nbr_cnt, ui *order,
                                                        TreeNode *& tree_node, ui curr_idx, VertexID data_vtx, size_t workload){

        //TreeNode curr_node = tree_node[curr_idx];

        if(tree_node[curr_idx].children_count_ == 0){
            return 1;
        }        
        

        for(ui i = 0; i < tree_node[curr_idx].children_count_; i++){
            
            VertexID child_id = tree_node[curr_idx].children_[i];

            ui* valid_nbrs = new ui[max_valid_nbr_cnt];
            ui nbr_cnt, valid_nbr_cnt = 0;
            
            ui query_vtx_degree = query_graph->getVertexDegree(child_id);

            const ui* nbrs = data_graph -> getVertexNeighbors(data_vtx, nbr_cnt);
            

            for(ui j = 0; j < nbr_cnt; j++){
                if(data_graph->getVertexDegree(nbrs[j]) >= query_vtx_degree){
                    valid_nbrs[valid_nbr_cnt++] = nbrs[j];
                }
            }

            //workload += valid_nbr_cnt;

            for(ui k = 0; k < valid_nbr_cnt; k++){
                workload *= calculateWorkLoad(data_graph, query_graph, max_valid_nbr_cnt, order, tree_node, child_id, valid_nbrs[k], workload);
            }

            delete[] valid_nbrs;

        }

        return workload;
}


size_t* LoadBalancer::workloadEstimator(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node){
        
        std::cout << "Estimating workload" << std::endl;

        size_t* est_workload_array = new size_t[candidates_count[0]];

        ui max_valid_nbr_cnt = data_graph -> getGraphMaxDegree();

        /*for(ui i = 0; i < query_graph->getVerticesCount(); i++){
            if(candidates_count[i] < max_valid_nbr_cnt){
                max_valid_nbr_cnt = candidates_count[i];
            }
        }*/

        for(ui i = 0; i < candidates_count[0]; i++){

            size_t workload = 0;
            
            est_workload_array[i] = calculateWorkLoad(data_graph, query_graph, max_valid_nbr_cnt, order, tree_node, 0, candidates[0][i], workload);

            // if(i % 1000 == 0){
            //     std::cout << "Done for : " << i << std::endl;
            // }
            
        }                                                   


        return est_workload_array;

}


void LoadBalancer::writeInCsvFile(int size, int rank, size_t* work_est_array, size_t* org_cnt_array, ui array_length, ui work_est_idx, ui org_cnt_idx){

    std::ofstream outfile;
    outfile.open ("comparison_result/P" + std::to_string(size) + "amazon_" + std::to_string(rank) + "_comparison.csv");
    outfile << "Estimated_Count,Real_Count," << std::endl;

    for(ui i = 0; i < array_length; i++){

        outfile << work_est_array[work_est_idx + i] << "," << org_cnt_array[org_cnt_idx + i] << "," << std::endl;

    }
    outfile.close();
}