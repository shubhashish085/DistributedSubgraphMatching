#include "ParallelEnumeration.h"
#include "Enumeration.h"
#include "GeneratingFilterPlan.h"
#include "wtime.h"
#include <omp.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#define NUM_THREADS 8
#define PAD 8


size_t** ParallelEnumeration::exploreWithDynamicLoadBalance(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                            TreeNode *& tree, size_t thread_output_limit_num, size_t &call_count, int &thread_count){

    std::cout << " ################## explore parallel ##################" << std::endl;

    size_t** embedding_cnt_array = new size_t * [thread_count];
    double* thread_wise_time = new double [thread_count];

    for(ui i = 0; i < thread_count; i++){
        embedding_cnt_array[i] = new size_t [PAD];
        embedding_cnt_array[i][0] = 0;
    }

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

    int max_depth = query_graph->getVerticesCount();
    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    VertexID start_vertex = order[0];

    omp_set_num_threads(thread_count);

    ui par_loop_idx = 0;


#pragma omp parallel
    {
        double start_time = wtime(), end_time;

        int th_id = omp_get_thread_num();
        int cur_depth = 0;

        VertexID **valid_candidate = new ui *[max_depth];
        ui *idx = new ui[max_depth];
        ui *idx_count = new ui[max_depth];
        ui *embedding = new ui[max_depth];
        bool *visited_vertices = new bool[data_graph->getVerticesCount()];

        VertexID* intersection_result = new VertexID[max_candidate_count];
        VertexID* intersection_order = new VertexID[max_depth];

        std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);

        for (ui i = 0; i < max_depth; ++i) {
            valid_candidate[i] = new VertexID[max_candidate_count];
        }


#pragma omp for schedule(dynamic, 5)
        for(par_loop_idx = 0; par_loop_idx < candidates_count[start_vertex]; par_loop_idx++){

            cur_depth = 0;
            idx[cur_depth] = 0;
            idx_count[cur_depth] = 1;
            valid_candidate[cur_depth][0] = candidates[start_vertex][par_loop_idx];

            while (true) {
                while (idx[cur_depth] < idx_count[cur_depth]) {
                    VertexID u = order[cur_depth];
                    VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
                    embedding[u] = v;
                    visited_vertices[v] = true;
                    idx[cur_depth] += 1;

                    if (cur_depth == max_depth - 1) {
                        embedding_cnt_array[th_id][0] += 1;
                        visited_vertices[v] = false;
                        //Enumerate::printMatch(embedding, query_graph->getVerticesCount());
                        /*if (embedding_cnt_array[th_id][0] >= thread_output_limit_num) {
                            goto EXIT;
                        }*/
                    } else {
                        call_count += 1;
                        cur_depth += 1;
                        idx[cur_depth] = 0;
                        Enumerate::generateValidCandidatesWithSetIntersectionByOrdering(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                                                        visited_vertices,tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);
                    }
                }

                cur_depth -= 1;
                if (cur_depth < 0)
                    break;
                else
                    visited_vertices[embedding[order[cur_depth]]] = false;
            }

            end_time = wtime();
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

        thread_wise_time[th_id] = end_time - start_time;

    }

    std::cout << "----------------------------------------------------------" << std::endl;

    std::cout << "Thread wise time spent for " << thread_count << " threads" << std::endl;

    for(ui i = 0; i < thread_count; i++){
        std::cout << "Thread ID : " << i << " Spent Time : " << thread_wise_time[i] << std::endl;
    }

    std::cout << "-----------------------------------------------------------" << std::endl;

    return embedding_cnt_array;
}