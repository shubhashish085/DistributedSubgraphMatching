
#ifndef GRAPH_H
#define GRAPH_H

#include <unordered_map>
#include <vector>
#include "types.h"

class Graph{

private:

    ui vertices_count;
    ui edges_count;
    ui labels_count;
    ui max_degree;
    ui max_label_frequency;

    ui* degrees;
    ui* offsets;
    VertexID * neighbors;
    ui* matching_order_idx;
    LabelID* labels;

    ui* reverse_index_offsets;
    ui* reverse_index;

    std::unordered_map<LabelID, ui> * neighborhood_label_count;
    std::unordered_map<LabelID, ui> labels_frequency;
    ui* labels_offsets;

public:

    Graph() {

        vertices_count = 0;
        edges_count = 0;
        labels_count = 0;
        max_degree = 0;
        max_label_frequency = 0;

        offsets = NULL;
        neighbors = NULL;
        labels = NULL;
        neighborhood_label_count = NULL;

        reverse_index_offsets = NULL;
        reverse_index = NULL;

        labels_offsets = NULL;
        labels_frequency.clear();

    }

    ~Graph() {
        delete[] offsets;
        delete[] neighbors;
        delete[] labels;
        delete[] reverse_index_offsets;
        delete[] reverse_index;
        delete[] neighborhood_label_count;
        delete[] labels_offsets;
    }

public:
    void loadGraphFromFile(const std::string& file_path);
    void loadGraphFromFileWithEdge(const std::string& file_path);
    void loadGraphFromFileWithoutStringConversion(const std::string& file_path);
    void loadGraphFromFileForWeakScaling(const std::string& file_path, ui division_factor);
    void loadGraphFromFileWithWeight(const std::string& file_path);
    void loadDirectedGraphFromFile(const std::string& file_path);
    void printGraphMetaData();
    void setMatchingOrderIndex(std::vector<ui> matching_order);

public:

    void BuildReverseIndex();

    const ui* getOffsets() const {
        return offsets;
    }

    const ui getLabelsCount() const {
        return labels_count;
    }

    const ui getVerticesCount() const {
        return vertices_count;
    }

    const ui getEdgesCount() const {
        return edges_count;
    }

    const ui getVertexDegree(const VertexID id) const {
        return offsets[id + 1] - offsets[id];
    }

    const ui getGraphMaxDegree() const {
        return max_degree;
    }

    const ui getGraphMaxLabelFrequency() const {
        return max_label_frequency;
    }

    const ui getLabelsFrequency(const LabelID label) const {
        return labels_frequency.find(label) == labels_frequency.end() ? 0 : labels_frequency.at(label);
    }

    const LabelID getVertexLabel(const VertexID id) const {
        return labels[id];
    }

    const ui * getVerticesByLabel(const LabelID id, ui& count) const {
        count = reverse_index_offsets[id + 1] - reverse_index_offsets[id];
        return reverse_index + reverse_index_offsets[id];
    }

    const ui * getNeighborsByLabel(const VertexID id, const LabelID label, ui& count) const {
        ui offset = id * labels_count + label;
        count = labels_offsets[offset + 1] - labels_offsets[offset];
        return neighbors + labels_offsets[offset];
    }

    const std::unordered_map<LabelID, ui>* getVertexNLF(const VertexID id) const {
        return neighborhood_label_count + id;
    }

    bool checkEdgeExistence(const VertexID u, const VertexID v, const LabelID u_label) const {
        ui count = 0;
        const VertexID* neighbors = getNeighborsByLabel(v, u_label, count);
        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    ui * getVertexNeighbors(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id]; // used for neighbor count
        return neighbors + offsets[id];
    }

    void getNeighborCount(const VertexID id, ui& count) const {
        count = offsets[id + 1] - offsets[id];
    }

    ui * getMatchingOrderIndex() const {
        return matching_order_idx;
    }

    std::unordered_map<LabelID, ui> * getNeighborhoodLabelCount(){
        return neighborhood_label_count;
    }

    const std::vector<VertexID> getVerticesOfLabel(LabelID labelId) const {

        std::vector<VertexID> candidate_vertices;

        for(VertexID i = 0; i < vertices_count; i++){
            if(labels[i] == labelId){
                candidate_vertices.push_back(i);
            }
        }

        return candidate_vertices;
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) {
            std::swap(u, v);
        }
        ui count = 0;
        const VertexID* neighbors =  getVertexNeighbors(v, count);

        int begin = 0;
        int end = count - 1;
        while (begin <= end) {
            int mid = begin + ((end - begin) >> 1);
            if (neighbors[mid] == u) {
                return true;
            }
            else if (neighbors[mid] > u)
                end = mid - 1;
            else
                begin = mid + 1;
        }

        return false;
    }

    void printGraphData(){

        std::cout << " >>>>> Degree <<<<<<" << std::endl;

        for(ui i = 0; i < vertices_count; i++){
            std::cout << degrees[i] << " ";
        }

        std::cout << std::endl;

        std::cout << " >>>>> Offset <<<<<<" << std::endl;

        for(ui i = 0; i < vertices_count + 1; i++){
            std::cout << offsets[i] << " ";
        }

        std::cout << std::endl;

        std::cout << " >>>>> Neighbors <<<<<<" << std::endl;

        for(ui i = 0; i < edges_count; i++){
            std::cout << neighbors[i] << " ";
        }

        std::cout << std::endl;

        std::cout << " >>>>> Labels <<<<<<" << std::endl;

        for(ui i = 0; i < vertices_count; i++){
            std::cout << labels[i] << " ";
        }

        std::cout << std::endl;
    }

};


#endif