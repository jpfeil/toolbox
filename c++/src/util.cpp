//
// Created by jpfeil on 11/14/19.
//
#include<string>
#include<iostream>
#include<algorithm>
#include<string_view>

#include<boost/regex.hpp>

std::string revcomp(std::string seq){
    auto lambda = [](const char c) {
        switch (c) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            case 'N':
                return 'N';
            default:
                throw std::domain_error("UTIL: Invalid nucleotide.");
        }
    };
    std::transform(seq.begin(), seq.end(), seq.begin(), lambda);
    std::reverse( seq.begin(), seq.end() );
    return seq;
}


std::string APOBEC(const std::string& kmer){
    boost::regex expr{"(TCA|TCT)"};
    auto lambda = [](boost::smatch const& match)->std::string{
        if (match.str(0) == "TCA"){
            return "TTA";
        } else if (match.str(0) == "TCT"){
            return "TTT";
        } else{
            return match.str(0);
        }
    };
    return boost::regex_replace(kmer, expr, lambda);
}


std::string APOBEC2(const std::string& kmer){
    boost::regex expr{"(TCA|TCT)"};
    auto lambda = [](boost::smatch const& match)->std::string{
        if (match.str(0) == "TCA"){
            return "TGA";
        } else if (match.str(0) == "TCT"){
            return "TGT";
        } else{
            return match.str(0);
        }
    };
    return boost::regex_replace(kmer, expr, lambda);
}

bool cas(int& value, int& old, int& nu){
    if (value != old){
        return false;
    }
    else{
        value = nu;
        return true;
    }
}
