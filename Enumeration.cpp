//
// Created by antu on 3/11/24.
//

#include "utilities.h"
#include "Enumeration.h"
#include <algorithm>
#include <map>
#include <fstream>

/* Tree list implementation */

void Enumerate::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                            ui *order, ui **candidates, ui *candidates_count) {

    //std::cout << " ################## generateValidCandidates ##################" << std::endl;

    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < tree[u].bn_count_; ++j) {
                VertexID u_nbr = tree[u].bn_[j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }

    std::cout << std::endl;

    std::cout << "Valid Candidates " << std::endl;
    for(ui i = 0; i <= depth; i++){
        for(ui j = 0; j < idx_count[i]; j++){
            std::cout << valid_candidate[i][j] << " " ;
        }
        std::cout << std::endl;
    }

}

/*
 * Generating Valid Candidates During Enumeration according to the previous vertex
 * */
void Enumerate::generateValidCandidatesWithCandidateCSR(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                    bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                    ui* candidate_offset, ui* candidate_csr){

    VertexID u = order[depth];
    ui neighbor_count = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> valid_candidate_map;
    std::map<ui, ui>::iterator search_result;

    for (ui i = 0; i < tree[u].bn_count_; i++){
        VertexID u_nbr = tree[u].bn_[i];
        VertexID u_nbr_v = embedding[u_nbr];

        VertexID* neighbors = data_graph ->getVertexNeighbors(u_nbr_v, neighbor_count);

        for(ui j = 0; j < neighbor_count; j++){
            VertexID v = neighbors[j];
            //isCandidateCheck
            for(ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++){
                if(candidate_csr[index] == u){

                    if(!visited_vertices[v] ){
                        search_result = valid_candidate_map.find(v);

                        if(search_result != valid_candidate_map.end()){
                            valid_candidate_map[v] = search_result -> second + 1;
                        }else{
                            valid_candidate_map[v] = 1;
                        }
                    }
                    break;
                }
            }
        }
    }

    std::map<ui, ui>::iterator it = valid_candidate_map.begin();

    ui bn_cout = tree[u].bn_count_;
    // Iterate through the map and print the elements
    while (it != valid_candidate_map.end()) {
        if(it -> second == bn_cout) {
            valid_candidate[depth][idx_count[depth]++] = it->first;
        }
        ++it;
    }

}

/*
 * Generating Valid Candidates During Enumeration according to the previous vertex using Set intersection and two pointer
 * */
void Enumerate::generateValidCandidatesWithSetIntersectionAndBinarySearch(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                      bool* visited_vertices,TreeNode *&tree, ui* order, ui* candidate_offset, ui* candidate_csr,
                                                                      VertexID* intersection_array, VertexID* intersection_order){
    VertexID u = order[depth];
    VertexID* neighbors;
    ui neighbor_count = 0;

    ui set_ints_length = 0, l_length = 0, max_intersection_boundary = INTMAX_MAX, min_intersection_boundary = 0, min_neighbor_count = INTMAX_MAX;
    int start_idx = 0, end_idx = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> intersection_map;
    std::map<ui, ui>::iterator search_result;

    ui bn_count = tree[u].bn_count_;

    for (ui i = 0; i < bn_count; i++){
        intersection_order[i] = embedding[tree[u].bn_[i]];
    }

    for (ui i = 0; i < bn_count; i++){
        neighbors = data_graph -> getVertexNeighbors(intersection_order[i], neighbor_count);

        if(neighbors[neighbor_count - 1] < max_intersection_boundary){
            max_intersection_boundary = neighbors[neighbor_count - 1];
        }

        if(neighbors[0] > min_intersection_boundary){
            min_intersection_boundary = neighbors[0];
        }

        if (neighbor_count < min_neighbor_count){
            ui temp = intersection_order[i];
            intersection_order[i] = intersection_order[0];
            intersection_order[0] = temp;
            min_neighbor_count = neighbor_count;
        }
    }


    for (ui i = 0; i < bn_count; i++){

        VertexID* neighbors = data_graph ->getVertexNeighbors(intersection_order[i], neighbor_count);

        if(neighbor_count > 15) {
            if (i == 0) {
                Utilities::binary_search_within_limit(neighbors, 0, neighbor_count, min_intersection_boundary,
                                                      max_intersection_boundary, start_idx, end_idx);
                std::copy(neighbors + start_idx, neighbors + end_idx, intersection_array);
                set_ints_length = end_idx - start_idx;
                l_length = set_ints_length;
            } else {
                Utilities::set_intersection_with_boundary_and_binary_search(intersection_array, l_length, neighbors,
                                                                            neighbor_count, set_ints_length,
                                                                            min_intersection_boundary,
                                                                            max_intersection_boundary);
                l_length = set_ints_length;
            }
        }else {
            if(i == 0){
                std::copy(neighbors, neighbors + neighbor_count, intersection_array);
                set_ints_length = neighbor_count;
                l_length = set_ints_length;
            }else{
                Utilities::set_intersection_tp(intersection_array, l_length, neighbors, neighbor_count, set_ints_length);
                l_length = set_ints_length;
            }
        }
    }

    for(ui i = 0; i < set_ints_length; i++){
        VertexID v = intersection_array[i];
        if(!visited_vertices[v]) {
            for (ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++) {
                if (candidate_csr[index] == u) {
                    valid_candidate[depth][idx_count[depth]++] = v;
                    break;
                }
            }
        }
    }

}

/*
 * Generating Valid Candidates During Enumeration according to the previous vertex using Set intersection and two pointer
 * */
void Enumerate::generateValidCandidatesWithSetIntersectionAndBoundary(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                     bool* visited_vertices,TreeNode *&tree, ui* order, ui* candidate_offset, ui* candidate_csr,
                                                                     VertexID* intersection_array, VertexID* intersection_order){
    VertexID u = order[depth];
    VertexID* neighbors;
    ui neighbor_count = 0;

    ui set_ints_length = 0, l_length = 0, max_intersection_boundary = INTMAX_MAX, min_intersection_boundary = 0, min_neighbor_count = INTMAX_MAX, start_idx = 0, end_idx = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> intersection_map;
    std::map<ui, ui>::iterator search_result;

    ui bn_count = tree[u].bn_count_;

    for (ui i = 0; i < bn_count; i++){
        intersection_order[i] = embedding[tree[u].bn_[i]];
    }

    for (ui i = 0; i < bn_count; i++){
        neighbors = data_graph -> getVertexNeighbors(intersection_order[i], neighbor_count);

        if(neighbors[neighbor_count - 1] < max_intersection_boundary){
            max_intersection_boundary = neighbors[neighbor_count - 1];
        }

        if(neighbors[0] > min_intersection_boundary){
            min_intersection_boundary = neighbors[0];
        }

        if (neighbor_count < min_neighbor_count){
            ui temp = intersection_order[i];
            intersection_order[i] = intersection_order[0];
            intersection_order[0] = temp;
            min_neighbor_count = neighbor_count;
        }
    }


    for (ui i = 0; i < bn_count; i++){

        VertexID* neighbors = data_graph ->getVertexNeighbors(intersection_order[i], neighbor_count);

        if(i == 0){
            start_idx = 0;
            end_idx = neighbor_count;

            while(start_idx < end_idx && neighbors[start_idx] < min_intersection_boundary){
                start_idx++;
            }

            while (start_idx < end_idx && neighbors[end_idx - 1] > max_intersection_boundary){
                end_idx--;
            }

            std::copy(neighbors + start_idx , neighbors + end_idx, intersection_array);
            set_ints_length = end_idx - start_idx;
            l_length = set_ints_length;
        }else{
            Utilities::set_intersection_with_boundary(intersection_array, l_length, neighbors, neighbor_count, set_ints_length, min_intersection_boundary, max_intersection_boundary);
            l_length = set_ints_length;
        }
    }

    for(ui i = 0; i < set_ints_length; i++){
        VertexID v = intersection_array[i];
        if(!visited_vertices[v]) {
            for (ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++) {
                if (candidate_csr[index] == u) {
                    valid_candidate[depth][idx_count[depth]++] = v;
                    break;
                }
            }
        }
    }

}


/*
 * Generating Valid Candidates During Enumeration according to the previous vertex using Set intersection and two pointer
 * */
void Enumerate::generateValidCandidatesWithSetIntersectionByOrdering(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                     bool* visited_vertices,TreeNode *&tree, ui* order, ui* candidate_offset, ui* candidate_csr,
                                                                     VertexID* intersection_array, VertexID* intersection_order){

    VertexID u = order[depth];
    ui neighbor_count = 0;

    ui set_ints_length = 0, l_length = 0, min_neighbor_count = INTMAX_MAX;

    idx_count[depth] = 0;

    std::map<ui, ui> intersection_map;
    std::map<ui, ui>::iterator search_result;

    ui bn_count = tree[u].bn_count_;

    for (ui i = 0; i < bn_count; i++){
        intersection_order[i] = embedding[tree[u].bn_[i]];
    }

    for (ui i = 0; i < bn_count; i++){
        data_graph ->getNeighborCount(intersection_order[i], neighbor_count);
        if (neighbor_count < min_neighbor_count){
            ui temp = intersection_order[i];
            intersection_order[i] = intersection_order[0];
            intersection_order[0] = temp;
            min_neighbor_count = neighbor_count;
        }
    }

    for (ui i = 0; i < bn_count; i++){

        VertexID* neighbors = data_graph ->getVertexNeighbors(intersection_order[i], neighbor_count);

        if(i == 0){
            std::copy(neighbors, neighbors + neighbor_count, intersection_array);
            set_ints_length = neighbor_count;
            l_length = set_ints_length;
        }else{
            Utilities::set_intersection_tp(intersection_array, l_length, neighbors, neighbor_count, set_ints_length);
            l_length = set_ints_length;
        }
    }

    for(ui i = 0; i < set_ints_length; i++){
        VertexID v = intersection_array[i];
        if(!visited_vertices[v]) {
            for (ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++) {
                if (candidate_csr[index] == u) {
                    valid_candidate[depth][idx_count[depth]++] = v;
                    break;
                }
            }
        }
    }

}



/*
 * Generating Valid Candidates During Enumeration according to the previous vertex using Set intersection and two pointer
 * */
void Enumerate::generateValidCandidatesWithSetIntersection_tp(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                           bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                           ui* candidate_offset, ui* candidate_csr, VertexID* intersection_array){

    VertexID u = order[depth];
    ui neighbor_count = 0;

    ui set_ints_length = 0, l_length = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> intersection_map;
    std::map<ui, ui>::iterator search_result;

    ui bn_count = tree[u].bn_count_;

    for (ui i = 0; i < tree[u].bn_count_; i++){
        VertexID u_nbr = tree[u].bn_[i];
        VertexID u_nbr_v = embedding[u_nbr];

        VertexID* neighbors = data_graph ->getVertexNeighbors(u_nbr_v, neighbor_count);

        if(i == 0){
            std::copy(neighbors, neighbors + neighbor_count, intersection_array);
            set_ints_length = neighbor_count;
            l_length = set_ints_length;
        }else{
            Utilities::set_intersection_tp(intersection_array, l_length, neighbors, neighbor_count, set_ints_length);
            l_length = set_ints_length;
        }
    }

    for(ui i = 0; i < set_ints_length; i++){
        VertexID v = intersection_array[i];
        if(!visited_vertices[v]) {
            for (ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++) {
                if (candidate_csr[index] == u) {
                    valid_candidate[depth][idx_count[depth]++] = v;
                    break;
                }
            }
        }
    }

}



/*
 * Generating Valid Candidates During Enumeration according to the previous vertex
 * */
void Enumerate::generateValidCandidatesWithSetIntersection(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                        ui* candidate_offset, ui* candidate_csr){

    VertexID u = order[depth];
    ui neighbor_count = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> intersection_map;
    std::map<ui, ui>::iterator search_result;

    ui bn_count = tree[u].bn_count_;

    for (ui i = 0; i < tree[u].bn_count_; i++){
        VertexID u_nbr = tree[u].bn_[i];
        VertexID u_nbr_v = embedding[u_nbr];

        VertexID* neighbors = data_graph ->getVertexNeighbors(u_nbr_v, neighbor_count);
        for(ui j = 0; j < neighbor_count; j++){
            VertexID v = neighbors[j];
            search_result = intersection_map.find(v);

            if(search_result != intersection_map.end()){
                intersection_map[v] = search_result -> second + 1;
            }else{
                intersection_map[v] = 1;
            }
        }
    }

    std::map<ui, ui>::iterator it = intersection_map.begin();

    while (it != intersection_map.end()){
        if(it -> second == bn_count){
            VertexID v = it->first;
            if(!visited_vertices[v]) {
                for (ui index = candidate_offset[v]; index < candidate_offset[v + 1]; index++) {
                    if (candidate_csr[index] == u) {
                        valid_candidate[depth][idx_count[depth]++] = it->first;
                        break;
                    }
                }
            }
        }
        ++it;
    }

}

/*
 * Generating Valid Candidates During Enumeration according to the previous vertex
 * */
void Enumerate::generateValidCandidatesWithBinarySearch(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                        ui* candidate_offset, ui* candidate_csr){

    VertexID u = order[depth];
    ui neighbor_count = 0;

    idx_count[depth] = 0;

    std::map<ui, ui> valid_candidate_map;
    std::map<ui, ui>::iterator search_result;

    for (ui i = 0; i < tree[u].bn_count_; i++){
        VertexID u_nbr = tree[u].bn_[i];
        VertexID u_nbr_v = embedding[u_nbr];

        VertexID* neighbors = data_graph ->getVertexNeighbors(u_nbr_v, neighbor_count);

        for(ui j = 0; j < neighbor_count; j++){
            VertexID v = neighbors[j];
            //isCandidateCheck
            int index = Utilities::binary_search(candidate_csr, candidate_offset[v], candidate_offset[v+1] - 1, u);
            if(index != -1) {
                if (!visited_vertices[v]) {
                    search_result = valid_candidate_map.find(v);

                    if (search_result != valid_candidate_map.end()) {
                        valid_candidate_map[v] = search_result->second + 1;
                    } else {
                        valid_candidate_map[v] = 1;
                    }
                }
            }
        }
    }

    std::map<ui, ui>::iterator it = valid_candidate_map.begin();

    ui bn_cout = tree[u].bn_count_;

    while (it != valid_candidate_map.end()) {
        if(it -> second == bn_cout) {
            valid_candidate[depth][idx_count[depth]++] = it->first;
        }
        ++it;
    }
}



void Enumerate::generateValidCandidatesForRecursive(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                        ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                        ui *order, ui **candidates, ui *candidates_count) {


    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i) {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v]) {
            bool valid = true;

            for (ui j = 0; j < tree[u].bn_count_; ++j) {
                VertexID u_nbr = tree[u].bn_[j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v)) {
                    valid = false;
                    break;
                }
            }

            if (valid) {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}


void Enumerate::printMatch(ui* embedding, ui max_depth){

    std::cout << "------------------ Matched -------------------" << std::endl;

    for(ui i = 0; i < max_depth; i++){
        std::cout << "(" << i << "," << embedding[i] << "),";
    }

    std::cout << std::endl;
}

void Enumerate::increment_vertex_participation(ui* embedding, ui max_depth, ui* vertex_participation){

    for(ui i = 0; i < max_depth; i++){
        vertex_participation[embedding[i]] += 1;
    }

}


void Enumerate::analyseAndWriteResult(const std::string& file_path, const Graph *data_graph, const Graph *query_graph,
                                      ui* vertex_participation_in_embedding){
    std::ofstream outputfile;
    outputfile.open(file_path);

    outputfile << "#VertexID    VertexDegree    VertexParticipationInEmbedding" << std::endl;

    for(ui i = 0; i < data_graph->getVerticesCount(); i++){
        outputfile << i << '\t' << data_graph->getVertexDegree(i) << '\t' << vertex_participation_in_embedding[i]<< std::endl;
        if(i % 1000){
            outputfile.flush();
        }
    }

    outputfile.flush();
    outputfile.close();
}


/*
 * Exploration and analysis of sequential algorithm
 */
size_t Enumerate::exploreAndAnalysis(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                          ui *candidates_count, ui *order, TreeNode *& tree, ui* vertex_participation_in_embedding,
                          size_t output_limit_num, size_t &call_count, const std::string& file_path) {

    std::cout << " ################## explore and analysis ##################" << std::endl;

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    call_count = 0;


    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    VertexID* intersection_order = new VertexID[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    std::fill(vertex_participation_in_embedding, vertex_participation_in_embedding + data_graph->getVerticesCount(), 0);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    std::cout << "Max candidate count : " << max_candidate_count << std::endl;
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    // candidate csr building
    ui* candidate_track = new ui[data_graph->getVerticesCount()];
    ui* candidate_offset = new ui[data_graph->getVerticesCount() + 1];
    ui candidate_csr_count = 0;

    std::fill(candidate_track, candidate_track + data_graph -> getVerticesCount(), 0);

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        candidate_csr_count += candidates_count[i];
    }

    ui* candidate_csr = new ui[candidate_csr_count];

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        for(ui j = 0; j < candidates_count[i]; j++){
            VertexID data_vertex = candidates[i][j];
            candidate_track[data_vertex]++;
        }
    }

    candidate_offset[0] = 0;

    for(ui i = 1; i < data_graph -> getVerticesCount() + 1; i++){
        candidate_offset[i] = candidate_offset[i - 1] + candidate_track[i - 1];
    }


    std::fill(candidate_track, candidate_track + data_graph -> getVerticesCount(), 0);

    for(ui i = 0; i < query_graph -> getVerticesCount(); i++){
        for(ui j = 0; j < candidates_count[i]; j++){
            VertexID data_vertex = candidates[i][j];
            candidate_csr[candidate_offset[data_vertex] + candidate_track[data_vertex]] = i;
            candidate_track[data_vertex]++;
        }
    }

    for (ui i = 0; i < data_graph->getVerticesCount() ; ++i) {
        std::sort(candidate_csr + candidate_offset[i], candidate_csr + candidate_offset[i + 1]); // sorting the query graph parent of every vertex
    }



    VertexID* intersection_result = new VertexID[max_candidate_count];
    ui intersection_length = 0;

    std::cout << "Candidate count of Start Vertex : " << candidates_count[start_vertex] << std::endl;
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    std::cout << "Entering Loop " << std::endl;

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            //std::cout << " v : " << v << std::endl;
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;
            call_count++;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                //printMatch(embedding, query_graph->getVerticesCount());
                //increment_vertex_participation(embedding, query_graph->getVerticesCount(), vertex_participation_in_embedding);
                if (embedding_cnt >= output_limit_num) {
                    std::cout << "Output Limit Exceeded" << std::endl;
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                /*generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                        visited_vertices, tree, order, candidates, candidates_count);*/
                /*generateValidCandidatesWithCandidateCSR(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, tree, order, candidates, candidates_count, candidate_offset, candidate_csr);*/
                /*generateValidCandidatesWithSetIntersection(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                        visited_vertices, tree, order, candidates, candidates_count, candidate_offset, candidate_csr);*/
                /*generateValidCandidatesWithSetIntersection_tp(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                           visited_vertices, tree, order, candidates, candidates_count, candidate_offset, candidate_csr, intersection_result);*/
                /*generateValidCandidatesWithSetIntersectionByOrdering(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                              visited_vertices,tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);*/
                /*generateValidCandidatesWithSetIntersectionAndBoundary(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                                     visited_vertices,tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);*/
                generateValidCandidatesWithSetIntersectionAndBinarySearch(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                                      visited_vertices,tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);

                /*generateValidCandidatesWithBinarySearch(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                        visited_vertices, tree, order, candidates, candidates_count, candidate_offset, candidate_csr);*/
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else if (cur_depth == 0){
            //std::cout << "Candidate covered : " << idx[cur_depth] << " of candidate count : " << idx_count[cur_depth] << " at level cur_depth : " << cur_depth << std::endl;
            visited_vertices[embedding[order[cur_depth]]] = false;
        }else {
            visited_vertices[embedding[order[cur_depth]]] = false;
        }
    }



    // Release the buffer.
    EXIT:
    //analyseAndWriteResult(file_path, data_graph, query_graph, vertex_participation_in_embedding);
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] valid_candidate[i];
    }
    delete[] vertex_participation_in_embedding;
    delete[] valid_candidate;
    delete[] intersection_result;



    std::cout << "Total Embedding Count : " << embedding_cnt << std::endl;
    std::cout << "Total Operation Count : " << call_count << std::endl;

    return embedding_cnt;
}


size_t Enumerate::explore(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order, TreeNode *& tree,
                                          size_t output_limit_num, size_t &call_count) {

    std::cout << " ################## explore ##################" << std::endl;

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    call_count = 0;


    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    std::cout << "Max candidate count : " << max_candidate_count << std::endl;
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }



    std::cout << "Candidate count of Start Vertex : " << candidates_count[start_vertex] << std::endl;
    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    std::cout << "Entering Loop " << std::endl;

    while (true) {
        //std::cout << " idx[cur_depth] : " << idx[cur_depth] << " idx_count[cur_depth] : " << idx_count[cur_depth] << std::endl;
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            //std::cout << " v : " << v << std::endl;
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;
            call_count++;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                //printMatch(embedding, query_graph->getVerticesCount());
                if (embedding_cnt >= output_limit_num) {
                    std::cout << "Output Limit Exceeded" << std::endl;
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, tree, order, candidates, candidates_count);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else if (cur_depth == 0){
            std::cout << "Candidate covered : " << idx[cur_depth] << " of candidate count : " << idx_count[cur_depth] << " at level cur_depth : " << cur_depth << std::endl;
            visited_vertices[embedding[order[cur_depth]]] = false;
        }else {
            visited_vertices[embedding[order[cur_depth]]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] valid_candidate[i];
    }

    delete[] valid_candidate;

    std::cout << " Total Embedding Count : " << embedding_cnt << std::endl;
    std::cout << "Total Operation Count : " << call_count << std::endl;

    return embedding_cnt;
}

void Enumerate::exploreRecursiveFashion(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *embedding,
                               ui curr_depth, ui max_depth, ui *order, ui* idx, ui* idx_count, bool* visited_vertices, VertexID **valid_candidate,
                               TreeNode *& tree, size_t &embedding_count, size_t &call_count){

    //std::cout << " Current Depth : " << curr_depth << std::endl;

    if(curr_depth == query_graph ->getVerticesCount() - 1){
        VertexID u = order[curr_depth];
        VertexID v = valid_candidate[curr_depth][idx[curr_depth]];
        embedding[u] = v;
        //printMatch(embedding, max_depth);
        embedding_count++;
        return;

    }else{

        VertexID u = order[curr_depth];
        VertexID v = valid_candidate[curr_depth][idx[curr_depth]];
        std::cout << " u: " << u << " v: " << v << std::endl;
        embedding[u] = v;
        visited_vertices[v] = true;
        idx[curr_depth]++;


        call_count += 1;
        curr_depth += 1;
        idx[curr_depth] = 0;
        generateValidCandidatesForRecursive(data_graph, curr_depth, embedding, idx_count, valid_candidate,
                                                visited_vertices, tree, order, candidates, candidates_count);

        //std::cout << " Valid Candidate Count : " << idx_count[curr_depth] << std::endl;

        for(ui i = 0; (curr_depth <= max_depth - 1) && (i < idx_count[curr_depth]); i++){
            exploreRecursiveFashion(data_graph, query_graph, candidates, candidates_count, embedding,
                                    curr_depth, max_depth, order, idx, idx_count, visited_vertices, valid_candidate,
                                    tree, embedding_count, call_count);
            idx[curr_depth]++;
        }

    }

    curr_depth -= 1;
    if(curr_depth == 0){
        visited_vertices[embedding[order[curr_depth]]] = false;
        return;
    }else {
        visited_vertices[embedding[order[curr_depth]]] = false;
    }

}

void Enumerate::exploreWithoutCandidate(const Graph *data_graph, const Graph *query_graph, ui *order, ui *embedding,
                                        ui curr_depth, bool* visited_vertices, TreeNode *& tree, size_t &embedding_count, size_t &call_count) {

    call_count++;

    if(curr_depth > query_graph->getVerticesCount() - 1){
        embedding_count++;
        //printMatch(embedding, query_graph->getVerticesCount());
        return;
    }

    VertexID u = order[curr_depth];

    //determining v_candidates and candidate_count
    VertexID * v_candidates = new VertexID [data_graph->getGraphMaxLabelFrequency()];
    ui candidate_count = 0;

    if(curr_depth == 0 && tree[u].bn_count_ == 0){

        for(ui i = 0; i < data_graph -> getVerticesCount(); i++){
            if((data_graph ->getVertexLabel(i) == query_graph ->getVertexLabel(u) && data_graph ->getVertexDegree(i) >= query_graph ->getVertexDegree(u))){
                v_candidates[candidate_count++] = i;
            }
        }
    }

    for(ui i = 0; i < tree[u].bn_count_; i++){
        VertexID prev_v = embedding[tree[u].bn_[i]];
        ui nbrs_count = 0;
        VertexID * prev_v_nbrs = data_graph ->getVertexNeighbors(prev_v, nbrs_count);

        ui valid_nbrs_count = 0;
        for(ui j = 0; j < nbrs_count; j++){
            if((data_graph ->getVertexLabel(prev_v_nbrs[j]) == query_graph ->getVertexLabel(u) && data_graph ->getVertexDegree(prev_v_nbrs[j]) >= query_graph ->getVertexDegree(u))){
                prev_v_nbrs[valid_nbrs_count++] = prev_v_nbrs[j];
            }
        }

        //set intersection
        if(i == 0){
            Utilities::set_copy(prev_v_nbrs, valid_nbrs_count, v_candidates, candidate_count);
        }else{
            Utilities::set_intersection(prev_v_nbrs, valid_nbrs_count, v_candidates, candidate_count);
        }

    }

    //std::cout << "Candidates of " << u << " : " ;
    for(ui i = 0; i < candidate_count; i++){
        std::cout << v_candidates[i] << ", ";
    }
    std::cout << std::endl;


    for(ui i = 0; i < candidate_count; i++){

        VertexID v = v_candidates[i];

        if(!visited_vertices[v] && (data_graph ->getVertexLabel(v) == query_graph ->getVertexLabel(u) &&
            data_graph ->getVertexDegree(v) >= query_graph ->getVertexDegree(u))){
            embedding[u] = v;
            visited_vertices[v] = true;
            curr_depth++;
            exploreWithoutCandidate(data_graph, query_graph, order, embedding,
                    curr_depth, visited_vertices, tree, embedding_count, call_count);
            curr_depth--;
            visited_vertices[v] = false;
        }
    }

    if(curr_depth <= 0){
        return;
    }

}



