//
// Created by jpfeil on 10/22/19.
//

#ifndef VACCINATE_PARSER_H
#define VACCINATE_PARSER_H

#include<string>
#include<vector>

using namespace std;

struct FASTA{
    std::string name;
    std::string seq;

    FASTA(std::string p_name,
          std::string p_seq)
            : name(std::move(p_name)),
              seq(std::move(p_seq)) {
        // No additional steps
    }

    bool operator==(const FASTA& other){
        return (name == other.name && seq == other.seq);
    }
};


struct FASTQ{
    std::string name;
    std::string seq;
    std::string comment;
    std::string quality;

    FASTQ(std::string& p_name,
            std::string& p_seq,
            std::string& p_comment,
            std::string& p_quality)
            : name(std::move(p_name)),
              seq(std::move(p_seq)),
              comment(std::move(p_comment)),
              quality(std::move(p_quality)){
        // No additional steps
    }
};

void pfasta(const std::string& path, std::vector<FASTA> &fastas);
void mmfasta(const std::string& path, std::vector<FASTA> &fastas);
void pfastq(const std::string& path, std::vector<FASTQ> &fastqs);
void mmfastq(const std::string& path, std::vector<FASTQ> &fastqs);

#endif //VACCINATE_PARSER_H
