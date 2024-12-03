#include "ParallelEnumeration.h"
#include "Enumeration.h"
#include "GeneratingFilterPlan.h"
#include "wtime.h"
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#define NUM_THREADS 8
#define PAD 8

size_t *ParallelEnumeration::exploreWithEvenDegreeDist(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                                                       TreeNode *&tree, size_t thread_output_limit_num, size_t &call_count)
{

    std::cout << " ################## explore parallel ##################" << std::endl;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    size_t *embedding_cnt_array = new size_t[world_size];
    double *thread_wise_time = new double[world_size];

    for (ui i = 0; i < world_size; i++)
    {
        embedding_cnt_array[i] = 0;
    }

    int max_depth = query_graph->getVerticesCount();
    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();

    VertexID start_vertex = order[0];

    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    VertexID *intersection_result = new VertexID[max_candidate_count];
    VertexID *intersection_order = new VertexID[max_depth];
    bool *visited_vertices = new bool[data_graph->getVerticesCount()];

    VertexID **valid_candidate = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i){
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    ui *candidate_track = new ui[data_graph->getVerticesCount()];
    ui *candidate_offset = new ui[data_graph->getVerticesCount() + 1];
    ui candidate_csr_count = 0;

    std::fill(candidate_track, candidate_track + data_graph->getVerticesCount(), 0);

    for (ui i = 0; i < query_graph->getVerticesCount(); i++)
    {
        candidate_csr_count += candidates_count[i];
    }

    ui *candidate_csr = new ui[candidate_csr_count];

    for (ui i = 0; i < query_graph->getVerticesCount(); i++)
    {
        for (ui j = 0; j < candidates_count[i]; j++)
        {
            VertexID data_vertex = candidates[i][j];
            candidate_track[data_vertex]++;
        }
    }

    candidate_offset[0] = 0;

    for (ui i = 1; i < data_graph->getVerticesCount() + 1; i++)
    {
        candidate_offset[i] = candidate_offset[i - 1] + candidate_track[i - 1];
    }

    std::fill(candidate_track, candidate_track + data_graph->getVerticesCount(), 0);

    for (ui i = 0; i < query_graph->getVerticesCount(); i++)
    {
        for (ui j = 0; j < candidates_count[i]; j++)
        {
            VertexID data_vertex = candidates[i][j];
            candidate_csr[candidate_offset[data_vertex] + candidate_track[data_vertex]] = i;
            candidate_track[data_vertex]++;
        }
    }

    for (ui i = 0; i < data_graph->getVerticesCount(); ++i)
    {
        std::sort(candidate_csr + candidate_offset[i], candidate_csr + candidate_offset[i + 1]); // sorting the query graph parent of every vertex
    }

    if (world_rank == 0)
    {

        ui *cand_degree_offset = new ui[candidates_count[start_vertex] + 1];
        cand_degree_offset[0] = 0;
        ui *candidate_limit = new ui[world_size];

        for (ui j = 1; j < candidates_count[start_vertex] + 1; j++)
        {
            cand_degree_offset[j] = cand_degree_offset[j - 1] + data_graph->getVertexDegree(candidates[start_vertex][j - 1]);
        }

        ui total_degree = cand_degree_offset[candidates_count[start_vertex]];
        ui avg_degree = total_degree / world_size;

        ui last_thread_degree_offset = 0;
        ui process_idx = 0;

        for (ui j = 1; j < candidates_count[start_vertex] + 1; j++)
        {
            if (process_idx == world_size - 1)
            {
                candidate_limit[process_idx++] = candidates_count[start_vertex];
                last_thread_degree_offset = cand_degree_offset[candidates_count[start_vertex]];
                break;
            }
            else if (j == candidates_count[start_vertex])
            {
                candidate_limit[process_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }
            else if (cand_degree_offset[j] - last_thread_degree_offset >= avg_degree)
            {
                candidate_limit[process_idx++] = j;
                last_thread_degree_offset = cand_degree_offset[j];
            }
        }

        for (int rank = 1; rank < world_size; ++rank)
        {
            MPI_Send(candidate_limit, world_size, MPI_UNSIGNED, rank, 1, MPI_COMM_WORLD);
        }

        double start_time = wtime();

        ui lower_idx = 0;

        std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
        

        int cur_depth = 0;

        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidate_limit[world_rank] - lower_idx;
        std::copy(candidates[start_vertex] + lower_idx, candidates[start_vertex] + candidate_limit[world_rank],
                  valid_candidate[cur_depth]);

        while (true)
        {
            while (idx[cur_depth] < idx_count[cur_depth])
            {
                VertexID u = order[cur_depth];
                VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
                embedding[u] = v;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1)
                {
                    embedding_cnt_array[world_rank] += 1;
                    visited_vertices[v] = false;
                    if (embedding_cnt_array[world_rank] >= thread_output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count += 1;
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    Enumerate::generateValidCandidatesWithSetIntersectionByOrdering(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                                                    visited_vertices, tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);
                }
            }

            cur_depth -= 1;
            if (cur_depth < 0)
                break;
            else
                visited_vertices[embedding[order[cur_depth]]] = false;
        }

        MPI_Status status;
        ui process_embedding_count = 0;
        for (ui rank = 1; rank < world_size; rank++)
        {
            MPI_Recv(&process_embedding_count, 1, MPI_UNSIGNED, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            embedding_cnt_array[rank] = process_embedding_count;
        }

        size_t final_embedding_count = 0;
        for (ui i = 0; i < world_size; i++)
        {
            final_embedding_count += embedding_cnt_array[i];
        }

        std::cout << "Total Embedding Count : " << final_embedding_count << std::endl;
    }
    else
    {

        MPI_Status status;
        ui *candidate_limit = new ui[world_size];
        MPI_Recv(candidate_limit, world_size, MPI_UNSIGNED, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        double start_time = wtime();

        ui lower_idx = candidate_limit[world_rank - 1];

        std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
        

        int cur_depth = 0;

        idx[cur_depth] = 0;
        idx_count[cur_depth] = candidate_limit[world_rank] - lower_idx;
        std::copy(candidates[start_vertex] + lower_idx, candidates[start_vertex] + candidate_limit[world_rank],
                  valid_candidate[cur_depth]);

        ui process_embedding_count = 0;

        while (true)
        {
            while (idx[cur_depth] < idx_count[cur_depth])
            {
                VertexID u = order[cur_depth];
                VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
                embedding[u] = v;
                visited_vertices[v] = true;
                idx[cur_depth] += 1;

                if (cur_depth == max_depth - 1)
                {
                    process_embedding_count += 1;
                    visited_vertices[v] = false;
                    if (process_embedding_count >= thread_output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count += 1;
                    cur_depth += 1;
                    idx[cur_depth] = 0;
                    Enumerate::generateValidCandidatesWithSetIntersectionByOrdering(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                                                                    visited_vertices, tree, order, candidate_offset, candidate_csr, intersection_result, intersection_order);
                }
            }

            cur_depth -= 1;
            if (cur_depth < 0)
                break;
            else
                visited_vertices[embedding[order[cur_depth]]] = false;
        }

        MPI_Send(&process_embedding_count, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
    }

// Release the buffer.
EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] valid_candidate[i];
    }

    delete[] valid_candidate;

    return embedding_cnt_array;
}
