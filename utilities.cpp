//
// Created by antu on 3/15/24.
//

#include "utilities.h"
#include <iostream>

void Utilities::set_copy(VertexID *l_array, ui &l_count, VertexID *valid_candidate, ui &candidate_count) {

    candidate_count = 0;

    for(ui i = 0; i < l_count; i++){
        valid_candidate[candidate_count++] = l_array[i];
    }

    return;
}

void Utilities::set_intersection_with_boundary_and_binary_search(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
                                               ui &set_ints_length, ui min, ui max) {
    int i = 0, j = 0, len = r_length;
    set_ints_length = 0;

    binary_search_within_limit(r_array, 0, len, min, max, j, len);

    while (i < l_length && j < len) {

        if (result[i] < r_array[j]) {
            i++;
        }
        else if (r_array[j] < result[i]) {
            j++;
        }
        else
        {
            result[set_ints_length++] = r_array[j];
            i++;
            j++;
        }
    }
}


void Utilities::set_intersection_with_boundary(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
                                                    ui &set_ints_length, ui min, ui max) {
    int i = 0, j = 0;
    set_ints_length = 0;

    while(j < r_length && r_array[j] < min){
        j++;
    }

    while (r_length > j && r_array[r_length - 1] > max){
        r_length--;
    }


    while (i < l_length && j < r_length) {

        if (result[i] < r_array[j]) {
            i++;
        }
        else if (r_array[j] < result[i]) {
            j++;
        }
        else
        {
            result[set_ints_length++] = r_array[j];
            i++;
            j++;
        }
    }
}

void Utilities::set_intersection_with_maximum_bound(VertexID *result, ui l_length, VertexID *r_array, ui r_length,
                                    ui &set_ints_length, ui max) {
    int i = 0, j = 0;
    set_ints_length = 0;

    while (i < l_length && j < r_length) {

        if (result[i] < r_array[j]) {
            i++;
            if(result[i] > max){
                break;
            }
        }
        else if (r_array[j] < result[i]) {
            j++;
            if(r_array[j] > max){
                break;
            }
        }
        else
        {
            result[set_ints_length++] = r_array[j];
            i++;
            j++;
        }
    }
}

void Utilities::set_intersection_tp(VertexID * result, ui l_length, VertexID * r_array, ui r_length, ui& set_ints_length) {

    int i = 0, j = 0;
    set_ints_length = 0;

    while (i < l_length && j < r_length) {
        if (result[i] < r_array[j])
            i++;
        else if (r_array[j] < result[i])
            j++;
        else
        {
            result[set_ints_length++] = r_array[j];
            i++;
            j++;
        }
    }

}

void Utilities::set_intersection(VertexID *r_array, ui &r_count,
                                 VertexID *valid_candidate, ui &candidate_count) {

    if(candidate_count == 0){
        return;
    }
    VertexID * l_array = new VertexID[candidate_count];
    ui  l_count = candidate_count;

    for(ui i = 0; i < l_count; i++){
        l_array[i] = valid_candidate[i];
    }

    candidate_count = 0;

    for(ui i = 0; i < l_count; i++){
        bool shouldBeInserted = false;
        for(ui j = 0; j < r_count; j++){
            if(l_array[i] == r_array[j]){
                shouldBeInserted = true;
                break;
            }
        }

        if(shouldBeInserted){
            valid_candidate[candidate_count++] = r_array[i];
        }
    }

    delete[] l_array;

    return;
}

int Utilities::binary_search (ui* array, ui low_idx, ui high_idx, ui element){

    while (low_idx <= high_idx) {
        int mid_idx = low_idx + (high_idx - low_idx) / 2;

        if (array[mid_idx] == element)
            return mid_idx;

        if (array[mid_idx] < element)
            low_idx = mid_idx + 1;

        else
            high_idx = mid_idx - 1;
    }

    return -1;
}

void Utilities::binary_search_within_limit(VertexID * arr, int low, int high, ui min, ui max, int& start_idx, int& end_idx)
{
    start_idx = -1, end_idx = -1;
    int mid = low;
    int start = low, end = high;

    while (low <= high) {
        mid = low + (high - low) / 2;

        if (arr[mid] == min){
            start_idx = mid;
            break;
        }

        if (arr[mid] < min){
            low = mid + 1;
        }

        else{
            high = mid - 1;
        }
    }

    if(start_idx == -1){
        if(arr[mid] > min && mid - 1 >= start){
            start_idx = mid - 1;
        }else{
            start_idx = mid;
        }
    }

    low = start;
    high = end;
    mid = low;

    while (low <= high) {
        mid = low + (high - low) / 2;


        if (arr[mid] == max){
            end_idx = mid;
            break;
        }

        if (arr[mid] < max){
            low = mid + 1;
        }

        else{
            high = mid - 1;
        }
    }

    if(end_idx == -1){
        if(arr[mid] < max && mid + 1 <= end){
            end_idx = mid + 1;
        }else{
            end_idx = mid;
        }
    }

    if(mid + 1 <= end){
        end_idx = mid + 1;
    }

    if(start_idx > end_idx){
        start_idx = start;
        end_idx = end;
    }

}

