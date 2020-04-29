//
// Created by jpfeil on 10/23/19.
//

#ifndef VACCINATE_HAMMING_H
#define VACCINATE_HAMMING_H

#include<string>
#include<set>

using namespace std;

int hamming_dist(const string& str1, const string& str2);

int min_ham(const string& str,
            const vector<string>& background,
            const int& minimum);

int min_ham(const string& str,
            const vector<string>& background,
            const int& minimum,
            const int& start,
            const int& stop);

#endif //VACCINATE_HAMMING_H
