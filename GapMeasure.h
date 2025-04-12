#ifndef DISTRIBUTEDSUBGRAPHMATCHING_GAPMEASURE_H
#define DISTRIBUTEDSUBGRAPHMATCHING_GAPMEASURE_H

#include "graph.h"

class GapMeasure {

public:
    static long long measure_the_epsilon_gap(const Graph *data_graph);
};


#endif //DISTRIBUTEDSUBGRAPHMATCHING_GAPMEASURE_H