

#include <algorithm>
#include <vector>
#include "types.h"
#include "Automorphism.h"

ui* Automorphism::convert_to_adj_mat(ui vertices_count, const ui* offset, const ui* neighbors){

    ui size = vertices_count;

    ui* adj_mat = new ui[size * size];
    std::fill(adj_mat, adj_mat + (size * size), 0);

    for(ui i = 0; i < vertices_count + 1; i++){
        for(ui j = offset[i]; j < offset[i + 1]; j++){
            adj_mat[i * size + neighbors[j]] = 1; 
        }
    }

    return adj_mat;
}

void Automorphism::get_automorphisms(std::vector<std::vector<ui>> &Aut, ui* adj_mat, ui size){
    
    ui p[size];
    Aut.clear();
    for(ui i = 0; i < size; ++i) p[i] = i;
    do{
        bool tag = true;
        for(ui i = 0; i < size; ++i) {
            for(ui j = 0; j < size; ++j)
                if( adj_mat[i * size + j] != adj_mat[p[i] * size +  p[j]]) {
                    tag = false;
                    break;
                }
            if( !tag ) break;
        }
        if(tag) {
            std::vector<ui> tmp;
            tmp.clear();
            for(ui i = 0; i < size; ++i) tmp.push_back(p[i]);
            Aut.push_back(tmp);
        }
    } while( std::next_permutation(p, p + size) );
}

void Automorphism::restriction_integration_with_scheduling(ui*& schedule, ui size, std::vector< std::pair<ui, ui>>& ordered_pairs, 
                                                    std::map<ui, std::vector<std::pair<ui, ui>>>& schedule_restriction_map){

    bool first_found, second_found;
    for(ui i = 0; i < ordered_pairs.size(); i++){

        for(ui j = 0; j < size; j++){
            first_found = false;
            second_found = false;

            for(ui k = 0; k <= j; k++){
                if(schedule[k] == ordered_pairs[i].first){
                    first_found = true;
                }else if(schedule[k] == ordered_pairs[i].second){
                    second_found = true;
                }
            }

            if(first_found && second_found){
                
                std::map<ui, std::vector<std::pair<ui, ui>>>::iterator it = schedule_restriction_map.find(j);

                if(it != schedule_restriction_map.end()){
                    (it->second).push_back(ordered_pairs[i]);
                }else{
                    std::vector<std::pair<ui, ui>> pair_vtr;
                    pair_vtr.push_back(ordered_pairs[i]);
                    schedule_restriction_map.insert({j, pair_vtr});
                }

                break;
            }  
        }
    }
}

void Automorphism:: aggressive_optimize(std::vector< std::pair<ui, ui>>& ordered_pairs, ui* adj_mat, ui size) { 
    std::vector< std::vector<ui>> Aut;
    get_automorphisms(Aut, adj_mat, size);

    std::vector< std::pair<ui,ui> > L;
    L.clear();

    for(ui v = 0; v < size; ++v) { // iterate all elements in schedule
        std::vector< std::vector<ui> > stabilized_aut;
        stabilized_aut.clear();

        for(ui i = 0; i < Aut.size(); ++i) {
            std::vector<ui>& x = Aut[i];
            if( x[v] == v) {
                stabilized_aut.push_back(x);
            }
            else {
                ui x1  = v, x2 = x[v];
                if( x1 > x2) {
                    ui tmp = x1;
                    x1 = x2;
                    x2 = tmp;
                }
                bool tag = true;
                std::pair<ui,ui> cur = std::make_pair(x1, x2);
                for(ui j = 0; j < L.size(); ++j)
                    if( L[j]  == cur) {
                        tag = false;
                        break;
                    }
                if(tag) {
                    L.push_back(cur);
                }
            }
        }

        Aut = stabilized_aut;
    }
    
    ordered_pairs.clear(); // In GraphZero paper, this vector's name is 'L'

    // pairing the optimized pairs which will reduce computation
    for(ui i = 0; i < L.size(); ++i) {
        bool tag = true;
        for(ui j = 0; j < ordered_pairs.size(); ++j)
            if( L[i].second == ordered_pairs[j].second) {
                tag = false;
                if( L[i].first > ordered_pairs[j].first) ordered_pairs[j].first = L[i].first;
                break;
            }
        if(tag) ordered_pairs.push_back(L[i]);
    }
}