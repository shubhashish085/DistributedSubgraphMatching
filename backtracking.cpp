//
// Created by antu on 2/8/24.
//

#include "backtracking.h"
#include <map>
#include <queue>

bool
AlgorithmStore :: DFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited){

    std::cout << "------------------------ In DFSTraversal Method ------------------------"  << std::endl;

    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    ui neighbor_count = 0;
    VertexID * neighbors =  graph->getVertexNeighbors(vertex_to_be_explored, neighbor_count);

    for(ui iterator = 0; iterator < neighbor_count; iterator++){
        if(!visited[neighbors[iterator]]){
            visited[neighbors[iterator]] = true;
            matching_order.push_back(neighbors[iterator]);
            DFSTraversal(graph, neighbors[iterator], non_tree_edges, matching_order, visited);
        } else{
            std::pair<VertexID, VertexID> nte = std::make_pair(vertex_to_be_explored, neighbors[iterator]);
            non_tree_edges.push_back(nte);
        }
    }

    return true;
}

void
AlgorithmStore::bfsTraversal(const Graph *graph, VertexID root_vertex, TreeNode *&tree, VertexID *&bfs_order) {

    //std::cout << "########################## BFS Traversal ############################" << std::endl;

    ui vertex_num = graph->getVerticesCount();

    std::queue<VertexID> bfs_queue;
    std::vector<bool> visited(vertex_num, false);

    tree = new TreeNode[vertex_num];
    for (ui i = 0; i < vertex_num; ++i) {
        tree[i].initialize(vertex_num);
    }
    bfs_order = new VertexID[vertex_num];

    ui visited_vertex_count = 0;
    bfs_queue.push(root_vertex);
    visited[root_vertex] = true;
    tree[root_vertex].level_ = 0;
    tree[root_vertex].id_ = root_vertex;

    while(!bfs_queue.empty()) {
        const VertexID u = bfs_queue.front();
        bfs_queue.pop();
        bfs_order[visited_vertex_count++] = u;

        ui u_nbrs_count;
        const VertexID* u_nbrs = graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui i = 0; i < u_nbrs_count; ++i) {
            VertexID u_nbr = u_nbrs[i];

            if (!visited[u_nbr]) {
                bfs_queue.push(u_nbr);
                visited[u_nbr] = true;
                tree[u_nbr].id_ = u_nbr;
                tree[u_nbr].parent_ = u;
                tree[u_nbr].level_ = tree[u] .level_ + 1;
                tree[u].children_[tree[u].children_count_++] = u_nbr;
            }
        }
    }

    //std::cout << "########################## End BFS Traversal ############################" << std::endl;

}


bool
AlgorithmStore :: BFSTraversal(const Graph *graph, VertexID start_node, std::vector<std::pair<VertexID , VertexID>>& non_tree_edges, std::vector<ui>& matching_order,
                               std::vector<bool>& visited, int* parent_vtr){

    std::cout << "------------------------ In BFSTraversal Method ------------------------"  << std::endl;
    visited[start_node] = true;
    VertexID vertex_to_be_explored = start_node;
    std::queue<VertexID> bfs_queue;
    bfs_queue.push(start_node);

    std::cout << "Start Node : " << start_node << std::endl;

    for(ui i = 0; i < graph -> getVerticesCount(); i++){
        parent_vtr[i] = -1;
    }

    std::vector<std::pair<VertexID , VertexID>> visited_edges;

    while(!bfs_queue.empty()){
        ui neighbor_count = 0;
        VertexID front_vertex = bfs_queue.front();
        bfs_queue.pop();
        VertexID * neighbors = graph ->getVertexNeighbors(front_vertex, neighbor_count);

        std::cout << "Neighbor count of " << front_vertex << " : " << neighbor_count << std::endl;

        matching_order.push_back(start_node);

        for(ui iterator = 0; iterator < neighbor_count; iterator++){

            if(!visited[neighbors[iterator]]){
                bfs_queue.push(neighbors[iterator]);
                matching_order.push_back(neighbors[iterator]);
                visited[neighbors[iterator]] = true;
                visited_edges.push_back(std::make_pair(front_vertex, neighbors[iterator]));
                parent_vtr[neighbors[iterator]] = start_node;
            }else{
                bool alreadyExists = false;
                for(ui i = 0; i < visited_edges.size(); i++){
                    if((visited_edges[i].first == front_vertex && visited_edges[i].second == neighbors[iterator]) ||
                    (visited_edges[i].first == neighbors[iterator] && visited_edges[i].second == front_vertex)){
                        alreadyExists = true;
                        break;
                    }
                }

                for(ui i = 0; i < non_tree_edges.size(); i++){
                    if((non_tree_edges[i].first == front_vertex && non_tree_edges[i].second == neighbors[iterator]) ||
                       (non_tree_edges[i].first == neighbors[iterator] && non_tree_edges[i].second == front_vertex)){
                        alreadyExists = true;
                        break;
                    }
                }


                if(!alreadyExists) {

                    std::pair <VertexID, VertexID> nte = std::make_pair(front_vertex, neighbors[iterator]);
                    non_tree_edges.push_back(nte);
                }
            }
        }
    }

    return true;
}


