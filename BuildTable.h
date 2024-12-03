//
// Created by kars1 on 5/23/24.
//

#ifndef DISTRIBUTEDSUBGRAPHMATCHING_BUILDTABLE_H
#define DISTRIBUTEDSUBGRAPHMATCHING_BUILDTABLE_H

#include "graph.h"

class BuildTable {

public:
    static size_t computeMemoryCostInBytes(const Graph *query_graph, const Graph *data_graph, ui *candidates_count);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_BUILDTABLE_H