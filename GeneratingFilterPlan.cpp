//
// Created by antu on 3/10/24.
//

#include "FilterVertices.h"
#include "GeneratingFilterPlan.h"
#include "backtracking.h"
#include <queue>
#include <algorithm>



VertexID GeneratingFilterPlan::selectCFLFilterStartVertex(const Graph *data_graph, const Graph *query_graph) {
    //min heap

    //std::cout << " #################### selectCFLFilterStartVertex ############## " << std::endl;
    auto rank_compare = [](std::pair<VertexID, double> l, std::pair<VertexID, double> r) {
        return l.second < r.second;
    };

    std::priority_queue<std::pair<VertexID, double>, std::vector<std::pair<VertexID, double>>, decltype(rank_compare)> rank_queue(rank_compare);

    // Compute the ranking.
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;

        //if (query_graph->get2CoreSize() == 0 || query_graph->getCoreValue(query_vertex) > 1) {
            LabelID label = query_graph->getVertexLabel(query_vertex);
            ui degree = query_graph->getVertexDegree(query_vertex);
            ui frequency = data_graph->getLabelsFrequency(label);
            double rank = frequency / (double) degree;
            rank_queue.push(std::make_pair(query_vertex, rank));
        //}
    }

    // Keep the top-3.
    while (rank_queue.size() > 3) {
        rank_queue.pop();
    }

    VertexID start_vertex = 0;
    double min_score = data_graph->getGraphMaxLabelFrequency() + 1;

    while (!rank_queue.empty()) {
        VertexID query_vertex = rank_queue.top().first;
        ui count;
        FilterVertices::computeCandidateWithNLF(data_graph, query_graph, query_vertex, count);
        double cur_score = count / (double) query_graph->getVertexDegree(query_vertex);

        if (cur_score < min_score) {
            start_vertex = query_vertex;
            min_score = cur_score;
        }
        rank_queue.pop();
    }

    return start_vertex;
}



void GeneratingFilterPlan::generateCFLFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                  VertexID *&order, int &level_count, ui *&level_offset) {

    //std::cout << " #################### GenerateCFLFilterPlan ############## " << std::endl;

    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID start_vertex = selectCFLFilterStartVertex(data_graph, query_graph);
    AlgorithmStore::bfsTraversal(query_graph, start_vertex, tree, order);

    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    level_count = -1;
    level_offset = new ui[query_vertices_num + 1];


    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        if (tree[u].level_ != level_count) {
            level_count += 1;
            level_offset[level_count] = 0;
        }

        level_offset[level_count] += 1;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];

            if (tree[u].level_ == tree[u_nbr].level_) {
                if (order_index[u_nbr] < order_index[u]) {
                    tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                }
                else {
                    tree[u].fn_[tree[u].fn_count_++] = u_nbr;
                }
            }
            else if (tree[u].level_ > tree[u_nbr].level_) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
            }
            else {
                tree[u].under_level_[tree[u].under_level_count_++] = u_nbr;
            }
        }
    }

    level_count += 1;

    ui prev_value = 0;
    for (ui i = 1; i <= level_count; ++i) {
        ui temp = level_offset[i];
        level_offset[i] = level_offset[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset[0] = 0;
}

ui* GeneratingFilterPlan::generateLoadBalacePlan(const Graph *data_graph, const Graph *query_graph, VertexID *&order, ui** candidates, ui* candidate_count, ui* candidate_offset, ui* candidate_csr) {

    ui* predicted_cardinality = new ui[data_graph -> getVerticesCount()];

    for(ui i = 0; i < data_graph -> getVerticesCount(); i++){
        predicted_cardinality[i] = 0;
    }

    VertexID start_vertex = order[0];
    VertexID second_vertex = order[1];
    ui neighbor_count = 0, count = 0;

    for (ui j = 0; j < candidate_count[start_vertex]; j++){

        count = 0;
        VertexID* neighbors = data_graph->getVertexNeighbors(candidates[start_vertex][j], neighbor_count);
        for(ui k = 0; k < neighbor_count; k++){
            VertexID neighbor = neighbors[k];

            for (ui i = candidate_offset[neighbor]; i < candidate_offset[neighbor + 1]; i++){
                if(candidate_csr[i] == second_vertex) {
                    count++;
                    break;
                }
            }
        }

        predicted_cardinality[candidates[start_vertex][j]] = count;

    }

    return predicted_cardinality;

}


void GeneratingFilterPlan::generateCFLFilterPlanForDirectedGraph(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                 VertexID *&order, int &level_count, ui *&level_offset) {

    //std::cout << " #################### GenerateCFLFilterPlanForDirectedGraph ############## " << std::endl;

    ui query_vertices_num = query_graph->getVerticesCount();
    VertexID start_vertex = selectCFLFilterStartVertex(data_graph, query_graph);
    AlgorithmStore::bfsTraversal(query_graph, start_vertex, tree, order);

    bool* visited = new bool [query_vertices_num];
    std::fill(visited, visited + query_vertices_num, false);

    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    level_count = -1;
    level_offset = new ui[query_vertices_num + 1];

    for (ui i = 0; i < query_vertices_num; ++i){

        VertexID u = order[i];

        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        if (tree[u].level_ != level_count) {
            level_count += 1;
            level_offset[level_count] = 0;
        }

        level_offset[level_count] += 1;

    }

    level_count += 1;

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID u = order[i];

        for(ui j = 0; j < i; j++){
            VertexID prev_u = order[j];

            ui prev_u_nbrs_count;
            const VertexID* prev_u_nbrs = query_graph->getVertexNeighbors(prev_u, prev_u_nbrs_count);

            for(ui k = 0; k < prev_u_nbrs_count; k++){
                if(prev_u_nbrs[k] == u){
                    if(tree[prev_u].level_ == tree[u].level_){
                        tree[prev_u].fn_[tree[prev_u].fn_count_++] = u;
                        tree[u].bn_[tree[u].bn_count_++] = prev_u;
                    }
                    if(tree[prev_u].level_ < tree[u].level_){
                        tree[prev_u].under_level_[tree[prev_u].under_level_count_++] = u;
                        tree[u].bn_[tree[u].bn_count_++] = prev_u;
                    }
                }
            }

        }


    }

    ui prev_value = 0;
    for (ui i = 1; i <= level_count; ++i) {
        ui temp = level_offset[i];
        level_offset[i] = level_offset[i - 1] + prev_value;
        prev_value = temp;
    }
    level_offset[0] = 0;

    std::cout << "Tree : " << std::endl;
    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        std::cout << "Back Neighbors of vertex " << i << " : " << std::endl;
        for(ui j = 0; j < tree[i].bn_count_; j++){
            std::cout << tree[i].bn_[j] << " " ;
        }
        std::cout << std::endl;

        std::cout << "Forward Neighbors of vertex " << i << " : " << std::endl;
        for(ui j = 0; j < tree[i].fn_count_; j++){
            std::cout << tree[i].fn_[j] << " " ;
        }
        std::cout << std::endl;

        std::cout << "Under Level of vertex " << i << " : " << std::endl;
        for(ui j = 0; j < tree[i].under_level_count_; j++){
            std::cout << tree[i].under_level_[j] << " " ;
        }
        std::cout << std::endl;
    }

}