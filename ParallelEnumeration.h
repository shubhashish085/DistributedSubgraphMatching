#ifndef DISTRIBUTEDSUBGRAPHMATCHING_PARALLELENUMERATION_H
#define DISTRIBUTEDSUBGRAPHMATCHING_PARALLELENUMERATION_H


#include "types.h"
#include "graph.h"


class ParallelEnumeration {

public:
    static size_t* exploreWithEvenDegreeDist(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                        TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count);
    

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_PARALLELENUMERATION_H