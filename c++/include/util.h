//
// Created by jpfeil on 11/14/19.
//

#ifndef VACCINATE_UTIL_H
#define VACCINATE_UTIL_H

#include<iostream>
#include<fstream>
#include<string>
#include<mutex>
#include<unordered_set>
#include<unordered_map>
#include<vector>

#include<tsl/sparse_map.h>
#include<tsl/sparse_set.h>

#include<vaccinaTE/parser.h>

using namespace std;

struct KmerDB {
    explicit KmerDB()= default;

    explicit KmerDB(size_t allocate){
        preallocate = allocate;
    };

    explicit KmerDB(std::vector<FASTA> input){
        fastas = std::move(input);
    }

    void add(std::string& kmer){
        int len = kmer.size();
        if( kmers.find(len) != kmers.end()){
            mutx.lock();
            kmers[len].insert(kmer);
            mutx.unlock();

        } else{
            mutx.lock();
            kmers[len] = {kmer};
            kmers[len].reserve(preallocate);
            mutx.unlock();

            mutx.lock();
            kmer_lengths.push_back(len);
            mutx.unlock();
        }
    }

    bool contains(const std::string& kmer){
        int len = kmer.size();
        if( kmers.find(len) != kmers.end()){
            return kmers[len].find(kmer) != kmers[len].end();
        } else{
            return false;
        }
    }

    unordered_set<std::string> get_db(int& len){
        if (kmers.find(len) != kmers.end()) {
            return kmers[len];
        } else{
            return unordered_set<std::string>{};
        }
    }

    vector<int> get_lengths(){
        return kmer_lengths;
    }

    void generate(int length){
        if( kmers.find(length) == kmers.end()) {
            kmer_lengths.push_back(length);
            kmers[length].reserve(fastas.size() * length);
        }
        for (const FASTA& fa: fastas) {
            if (fa.seq.size() > length) {
                for (size_t i = 0; i < fa.seq.size() - length + 1; ++i) {
                    std::string kmer = fa.seq.substr(i, length);
                    add(kmer);
                }
            }
        }
    }

private:
    unordered_map<int, unordered_set<std::string>> kmers{};
    vector<int> kmer_lengths{};
    size_t preallocate = 100000;
    std::mutex mutx;
    vector<FASTA> fastas;
};


struct ProbeDB {
    tsl::sparse_map<std::string, std::string> probes{};
    tsl::sparse_set<std::string> multis;

    explicit ProbeDB()= default;

    explicit ProbeDB(std::vector<FASTA>& input){
        fastas = std::move(input);
    }

    void add(const std::string& probe, const std::string& name){
        if( multis.find(probe) != multis.end() ){
            // Do nothing

        } else if (probes.find(probe) == probes.end()){
            // If probe is not in probes map, then
            // add it. It may be a unique kmer

            probes[probe] = name;

        } else{
            // Probe was already inserted, so remove
            // it because it is not unique

            probes.erase(probe);
            multis.insert(probe);
        }
    }

    void generate(int length){
        std::string tmp{};
        for (const FASTA& fa: fastas) {
            if (fa.seq.size() > length) {
                for (size_t i = 0; i < fa.seq.size() - length + 1; ++i) {
                    tmp = fa.seq.substr(i, length);
                    add(tmp, fa.name);
                }
            }
        }
    }

    void write(std::string& path){
        ofstream myfile(path);
        if (myfile.is_open()) {
            for( const auto& pair: probes){
                myfile << pair.first << "\t" << pair.second << "\n";
            }
        }
        myfile.close();
    }

private:
    vector<FASTA> fastas;
};

std::string revcomp(std::string seq);

int cas(int& value, int& old, int& nu);

#endif //VACCINATE_UTIL_H
