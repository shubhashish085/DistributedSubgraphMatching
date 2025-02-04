#ifndef DISTRIBUTEDSUBGRAPHMATCHING_AUTOMORPHISM_H
#define DISTRIBUTEDSUBGRAPHMATCHING_AUTOMORPHISM_H


#include "graph.h"
#include <vector>
#include <map>


class Automorphism {

public:
    static ui* convert_to_adj_mat(ui vertices_count, const ui* offset, const ui* neighbors);
    static void get_automorphisms(std::vector<std::vector<ui>> &Aut, ui* adj_mat, ui size);
    static void aggressive_optimize(std::vector< std::pair<ui, ui>>& ordered_pairs, ui* adj_mat, ui size);
    static void restriction_integration_with_scheduling(ui*& schedule, ui size, std::vector< std::pair<ui, ui>>& ordered_pairs, 
                                                    std::map<ui, std::vector<std::pair<ui, ui>>>& schedule_restriction_map);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_AUTOMORPHISM_H
