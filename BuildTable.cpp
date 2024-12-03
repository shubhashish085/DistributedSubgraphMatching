//
// Created by kars1 on 5/23/24.
//

#include "BuildTable.h"

size_t BuildTable::computeMemoryCostInBytes(const Graph *query_graph, const Graph *data_graph, ui *candidates_count) {
    size_t memory_cost_in_bytes = 0;
    size_t per_element_size = sizeof(ui);

    ui query_vertices_num = query_graph->getVerticesCount();
    for (ui i = 0; i < query_vertices_num; ++i) {
        memory_cost_in_bytes += candidates_count[i] * per_element_size;
    }

    //adding memory of candidate csr
    memory_cost_in_bytes = memory_cost_in_bytes * 2;

    //adding memory of candidate offset
    memory_cost_in_bytes += data_graph->getVerticesCount() * per_element_size;

    return memory_cost_in_bytes;
}