#ifndef DISTRIBUTEDSUBGRAPHMATCHING_UTILITIES_H
#define DISTRIBUTEDSUBGRAPHMATCHING_UTILITIES_H

#include "types.h"

class Utilities {

public:
    static void set_copy (VertexID* l_array, ui &l_count, VertexID* valid_candidate, ui &candidate_count);

    static void set_intersection_tp(VertexID * result, ui l_length, VertexID * r_array, ui r_length, ui& set_ints_length);

    static void set_intersection_with_maximum_bound(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
            ui &set_ints_length, ui max);

    static void set_intersection_with_boundary(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
            ui &set_ints_length, ui min, ui max);

    static void set_intersection_with_boundary_and_binary_search(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
                                                                 ui &set_ints_length, ui min, ui max);

    static void set_intersection (VertexID* r_array, ui &r_count, VertexID* valid_candidate, ui &candidate_count);

    static int binary_search (ui* array, ui low_idx, ui high_idx, ui element);

    static void binary_search_within_limit(VertexID * arr, int low, int high, ui min, ui max, int& start_idx, int& end_idx);

};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_UTILITIES_H