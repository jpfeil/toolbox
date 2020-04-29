//
// Created by jpfeil on 11/14/19.
//

#ifndef VACCINATE_KMER_H
#define VACCINATE_KMER_H

#include<mutex>
#include<string>
#include<vector>
#include<unordered_set>
#include<unordered_map>
#include<string_view>

#include<vaccinaTE/parser.h>
#include<vaccinaTE/util.h>

struct Counter {
    std::unordered_map<std::string, int> count{};
    std::vector<size_t> kmer_lengths{};

    explicit Counter(const std::vector<FASTA>& kmers): kmers{kmers} {
        for (const auto& fasta: kmers){
            auto it = std::find(kmer_lengths.begin(),
                                kmer_lengths.end(),
                                fasta.seq.size());
            if (it == kmer_lengths.end()) {
                kmer_lengths.push_back(fasta.seq.size());
            }
            count[fasta.seq] = 0;
        }
    }

    void increment(std::string& key){
        if(count.find(key) != count.end()){
            bool done = false;
            while( !done ){
                int old = count[key];
                int nu = old + 1;
                if ( cas(count[key], old, nu) ){
                    done = true;
                }
            }
        }
        else{
            count[key] = 1;
        }
    }

private:
    const vector<FASTA>& kmers;
};


void populate(const int& length,
              const string& sequence,
              unordered_set<string>& container);


void kmer_map(const string& name,
              const int& length,
              const string& sequence,
              unordered_map<string, vector<string>>& kmap);


void kmer_count(const int& length,
                const string& sequence,
                unordered_map<string, int>& counter);


void kmer_finder(const std::string& read,
                 const std::string& quality,
                 const std::string& readtype,
                 const std::string& libtype,
                 const int &min_quality,
                 Counter &counter);


void count_reads(const string &path,
                 const string &libtype,
                 const string &readtype,
                 Counter& counter,
                 int &total,
                 const int &min_quality);

void probe_reads(const string &path,
                 const string &libtype,
                 const string &readtype,
                 Counter& vcounter,
                 Counter& mcounter,
                 Counter& ncounter,
                 int &total,
                 const int& min_quality);

unordered_set<string> batch_populate(const int& length,
                           const vector<FASTA>& sequences);


#endif //VACCINATE_KMER_H
