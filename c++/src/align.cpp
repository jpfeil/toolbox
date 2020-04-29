//
// Created by jpfeil on 2/11/20.
//
#include<algorithm>
#include<iostream>
#include<cstdio>
#include<boost/multi_array.hpp>

typedef boost::multi_array<int, 2> matrix;

int gap_pen = -5;

int s(char a, char b){
    if(a == b){
        return 1;
    }
    else{
        return -5;
    }
}

void init(boost::multi_array<int, 2> H){
    // Set first column to zero
    for(size_t i=0; i < H.shape()[0]; ++i){
        H[i][0] = 0;
    }

    // Set first row to zero
    for(size_t j=0; j < H.shape()[1]; ++j){
        H[0][j] = 0;
    }
}


std::pair<int, int> fill(std::string& A, std::string& B, boost::multi_array<int, 2>& H) {
    int match, del, ins, max;
    max = 0;
    std::pair<int, int> index;
    // Fill in scoring matrix
    for (size_t i = 1; i <= A.size(); ++i) {
        for (size_t j = 1; j <= B.size(); ++j) {
            match = H[i - 1][j - 1] + s(A[i-1], B[j-1]);

            del = H[i - 1][j] + gap_pen;
            ins = H[i][j - 1] + gap_pen;

            std::vector<int> params = {match, ins, del, 0};
            H[i][j] = *std::max_element(params.begin(), params.end());

            if( H[i][j] > max ){
                index = std::make_pair<int, int>(i, j) ;
                max = H[i][j];
            }
        }
    }
    return index;
}


void print(boost::multi_array<int, 2> H) {
    FILE *fout = fopen("/tmp/align.tsv", "w");
    for (size_t i = 0; i < H.shape()[0]; ++i) {
        for (size_t j = 0; j < H.shape()[1]; ++j) {
            fprintf(fout, "%d\t", H[i][j]);
        }
        fprintf(fout, "\n");
    }
}


int sw(std::string& A, std::string& B){
    matrix H(boost::extents[ A.size() + 1 ][ B.size() + 1 ]);
    init(H);
    std::pair<int, int> max = fill(A, B, H);
    return H[max.first][max.second];
}
