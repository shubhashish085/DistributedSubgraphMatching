//
// Created by antu on 3/11/24.
//

#ifndef DISTRIBUTEDSUBGRAPHMATCHING_ENUMERATE_H
#define DISTRIBUTEDSUBGRAPHMATCHING_ENUMERATE_H

#include "graph.h"
#include "types.h"

class Enumerate {

public:
    static void printMatch(ui* embedding, ui max_depth);

    static void increment_vertex_participation(ui* embedding, ui max_depth, ui* vertex_participation);

    static void analyseAndWriteResult(const std::string& file_path, const Graph *data_graph, const Graph *query_graph,
                          ui* vertex_participation_in_embedding);

    static void exploreWithoutCandidate(const Graph *data_graph, const Graph *query_graph, ui *order, ui *embedding,
                            ui curr_depth, bool* visited_vertices, TreeNode *& tree, size_t &embedding_count, size_t &call_count);

    static void exploreRecursiveFashion(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *embedding,
                            ui curr_depth, ui max_depth, ui *order, ui* idx, ui* idx_count, bool* visited_vertices, VertexID **valid_candidate,
                            TreeNode *& tree, size_t &embedding_count, size_t &call_count);

    static size_t explore(const Graph *data_graph, const Graph *query_graph, ui **candidates, ui *candidates_count, ui *order,
                          TreeNode *& tree, size_t output_limit_num, size_t &call_count);

    static size_t exploreAndAnalysis(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                   ui *candidates_count, ui *order, TreeNode *& tree, ui* vertex_participation_in_embedding,
                                   size_t output_limit_num, size_t &call_count, const std::string& file_path);

    static void generateValidCandidates(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count);

    static void generateValidCandidatesWithCandidateCSR(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                   bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                        ui* candidate_offset, ui* candidate_csr);

    static void generateValidCandidatesWithBinarySearch(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                        bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                        ui* candidate_offset, ui* candidate_csr);

    static void generateValidCandidatesWithSetIntersection(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                           bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                           ui* candidate_offset, ui* candidate_csr);

    static void generateValidCandidatesWithSetIntersection_tp(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                              bool* visited_vertices, TreeNode *&tree, ui* order, ui **candidates, ui* candidates_count,
                                                              ui* candidate_offset, ui* candidate_csr, VertexID* intersection_array);

    static void generateValidCandidatesWithSetIntersectionByOrdering(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                     bool* visited_vertices, TreeNode *&tree, ui* order,
                                                              ui* candidate_offset, ui* candidate_csr, VertexID* intersection_array, VertexID* intersection_order);

    static void generateValidCandidatesWithSetIntersectionAndBoundary(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                      bool* visited_vertices,TreeNode *&tree, ui* order, ui* candidate_offset, ui* candidate_csr,
                                                                      VertexID* intersection_array, VertexID* intersection_order);

    static void generateValidCandidatesWithSetIntersectionAndBinarySearch(const Graph* data_graph, ui depth, ui* embedding, ui* idx_count, ui** valid_candidate,
                                                                      bool* visited_vertices,TreeNode *&tree, ui* order, ui* candidate_offset, ui* candidate_csr,
                                                                      VertexID* intersection_array, VertexID* intersection_order);

    static void generateValidCandidatesForRecursive(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                        ui **valid_candidate, bool *visited_vertices, TreeNode *&tree,
                                        ui *order, ui **candidates, ui *candidates_count);


};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_ENUMERATE_H
