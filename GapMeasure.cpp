#include <iostream>
#include <cmath>
#include "GapMeasure.h"


long long GapMeasure::measure_the_epsilon_gap(const Graph *data_graph){

    long long gap_distance = 0;

    for(ui i = 0; i < data_graph->getVerticesCount(); i++){
        for(ui j = data_graph->offsets[i]; j < data_graph->offsets[i + 1]; j++){
            if(i < data_graph-> neighbors[j]){
                gap_distance += (long long)std::abs((int)(data_graph->neighbors[j] - i));
            }
        }
    }

    std::cout << "The gap distance : " << gap_distance << std::endl;


    return gap_distance;
}


long long GapMeasure::measure_the_beta_gap(const Graph *data_graph){

    long long gap_distance = 0;
    VertexID max_value = 0;

    for(ui i = 0; i < data_graph->getVerticesCount(); i++){
        max_value = 0;        
        for(ui j = data_graph->offsets[i]; j < data_graph->offsets[i + 1]; j++){
            if(i < data_graph-> neighbors[j]){
                max_value = (VertexID)std::max((int)max_value, std::abs((int)(data_graph->neighbors[j] - i)));
            }
        }
        gap_distance += max_value;
    }

    std::cout << "The gap distance : " << gap_distance << std::endl;

    return gap_distance;
}

