//
// Created by jpfeil on 10/23/19.
//
#include<iostream>
#include<fstream>
#include<string>
#include<boost/program_options.hpp>
#include"parser.h"
#include"hamming.h"


void search(std::vector<FASTA> &epitopes,
            std::vector<std::string> &reference,
            std::map<std::string,int> &minham){

    for (FASTA epitope: epitopes){
        int ham = min_ham(epitope.seq, reference);
        minham[epitope.name] = ham;
    }
}


int main(int argc, char** argv){
    using namespace boost::program_options;
    bool is_help{};

    options_description description{"protham -r <reference> <Epitope FASTA>"};
    description.add_options()
            ("help,h", bool_switch(&is_help), "display a help dialog")
            ("reference,r", value<std::string>(), "Path to epitope peptide FASTA (Protein-space)");

    positional_options_description positional;
    positional.add("epitopes", 1);

    command_line_parser parser{ argc, argv };
    parser.options(description);
    parser.positional(positional);

    variables_map vm;
    try{
        auto parsed_result = parser.run();
        store(parsed_result, vm);
        notify(vm);
    } catch (const std::exception& e){
        std::cerr << e.what() << "\n";
        return -1;
    }

    if (is_help){
        std::cout << description;
        return 0;
    }

    if (vm["epitopes"].empty()){
        std::cerr << "You must provide an epitope FASTA.\n";
        return -1;
    }

    if (vm["reference"].empty()){
        std::cerr << "You must provide reference FASTA.\n";
        return -1;
    }

    const auto &path = vm["epitopes"].as<std::string>();
    std::vector<FASTA> epitopes{};
    const char* pointer = path.c_str();
    pfasta(pointer, epitopes);

    const auto &path2 = vm["reference"].as<std::string>();
    std::vector<FASTA> proteome{};
    const char* pointer2 = path2.c_str();
    pfasta(pointer2, proteome);

    std::vector<std::string> reference{};
    for (FASTA element: proteome){
        reference.push_back(element.seq);
    }

    std::map<std::string, int> minham{};
    search(epitopes, reference, minham);

    ofstream myfile("vaccinaTE-minham.tsv");
    if (myfile.is_open()){
        for (std::pair<std::string, int> element: minham){
            std::string name = element.first;
            int count = element.second;
            myfile << name << "\t" << count << "\n";
        }
    }
    myfile.close();
}

