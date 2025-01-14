#ifndef DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H
#define DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H


#include "types.h"
#include "graph.h"


class LoadBalancer {

public:
    static size_t* workloadEstimator(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node);
    static size_t calculateWorkLoad(const Graph *data_graph, const Graph *query_graph, ui max_valid_nbr_cnt, ui *order,
                                                        TreeNode *& tree_node, ui curr_idx, VertexID data_vtx, size_t workload);
    static void writeInCsvFile(int size, int rank, size_t* work_est_array, size_t* org_cnt_array, ui array_length, ui work_est_idx, ui org_cnt_idx);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H