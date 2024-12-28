#ifndef DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H
#define DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H


#include "types.h"
#include "graph.h"


class LoadBalancer {

public:
    static size_t* workloadEstimator(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node);
    static size_t calculateWorkLoad(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree_node, ui curr_idx, size_t& workload);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_LOADBALANCER_H