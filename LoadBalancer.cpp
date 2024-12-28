
#include "types.h"
#include "LoadBalancer.h"


size_t LoadBalancer::calculateWorkLoad(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node, ui curr_idx, size_t& workload){

        TreeNode curr_node = tree_node[curr_idx];

        if(curr_node.children_count_ == 0){

            return workload;
        }


        for(ui i = 0; i < curr_node.children_count_; i++){
            
            VertexID child_id = curr_node.children_[i];
            LabelID label_id = query_graph -> getVertexLabel(child_id);

            ui label_nbr_cnt;
            const ui* label_nbr = data_graph -> getNeighborsByLabel(child_id, label_id, label_nbr_cnt);

            workload += label_nbr_cnt;

            workload += calculateWorkLoad(data_graph, query_graph, candidates, candidates_count, order, tree_node, child_id, workload);

        }



}


size_t* LoadBalancer::workloadEstimator(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node){
        
        size_t* est_workload_array = new size_t[candidates_count[0]];

        for(ui i = 0; i < candidates_count[0]; i++){

            size_t workload = 0;
            
            est_workload_array[i] = calculateWorkLoad(data_graph, query_graph, candidates, candidates_count, order, tree_node, 0, workload); 
        }                                                   


        return est_workload_array;

}