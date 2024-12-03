//
// Created by antu on 3/10/24.
//

#ifndef DISTRIBUTEDSUBGRAPHMATCHING_FILTERVERTICES_H
#define DISTRIBUTEDSUBGRAPHMATCHING_FILTERVERTICES_H


#include "graph.h"
#include <vector>


class FilterVertices {

public:
    static bool LDFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count);
    static bool NLFFilter(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);
    static bool NLFAndLDFFilter(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);
    static bool CFLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order, TreeNode *&tree);

    static void computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                        ui &count, ui *buffer = NULL);

    static void computeCandidateWithLDF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                        ui &count, ui *buffer = NULL);

    static void generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                   VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                   ui *candidates_count, ui *flag, ui *updated_flag);

    static void pruneCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                ui *candidates_count, ui *flag, ui *updated_flag);
private:
    static void allocateBuffer(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);

    static void compactCandidates(ui** &candidates, ui* &candidates_count, ui query_vertex_num);

    static bool isCandidateSetValid(ui** &candidates, ui* &candidates_count, ui query_vertex_num);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_FILTERVERTICES_H
