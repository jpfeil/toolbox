#include<bits/stdc++.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>

#include<boost/iostreams/device/mapped_file.hpp>
#include<boost/iostreams/filtering_streambuf.hpp>
#include<boost/iostreams/filter/gzip.hpp>
#include<boost/interprocess/file_mapping.hpp>
#include<boost/iostreams/stream.hpp>

#include<vaccinaTE/parser.h>

using namespace std;
namespace bip = boost::interprocess;
namespace bio = boost::iostreams;


void pfasta(const string& path, vector<FASTA>& fastas){
    /*
     * Parses FASTA file and loads vector with FASTA objects
     */

    ifstream file(path,
                  ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    // Add gzip decompressor if needed
    if (path.substr(path.find_last_of('.') + 1) == "gz"){
        in.push(boost::iostreams::gzip_decompressor());
    }
    in.push(file);
    std::istream input{&in};
	if (!input.good()){
		std::cerr << "Error opening '" << path << "'. Bailing out." << std::endl;
	}
	std::string line{}, name{}, content{};
	while( std::getline( input, line ) ){
	    if ( line.empty() ){
	        cerr << "WARNING: corrupt line!" << endl;
	        name.clear();
	        content.clear();
	        continue;
	    }
	    if( line[0] == '>' ){
	        // Output complete entry
	        if( !name.empty() ) {
                fastas.emplace_back( FASTA{name, content} );
            }
	        // Save new entries header line
	        if ( !line.empty() ) {
	            name = line.substr(1);
            }
	        // Clear content for new entry
	        content.clear();
	        // Found sequence for previous name
	    } else if ( !name.empty() ){
	        content += line;
	    }
	}
	// Reached the end, so save final entry
	if ( !name.empty() ){
        fastas.emplace_back( name, content );
    }
}


void mmfasta(const string& path, vector<FASTA> &fastas){
    /*
     * Parses FASTA file and loads vector with FASTA objects
     */
    bio::stream<bio::mapped_file_source> file;
    file.open(bio::mapped_file_source(path));

    bio::filtering_streambuf<bio::input> in;

    // Add gzip decompressor if needed
    if (path.substr(path.find_last_of('.') + 1) == "gz"){
        in.push(bio::gzip_decompressor());
    }

    in.push(file);
    std::istream input{&in};
    if (!input.good()){
        std::cerr << "Error opening '" << path << "'. Bailing out." << std::endl;
    }

    std::string line{}, name{}, content{};
    while( std::getline( input, line ) ){
        if ( line.empty() ){
            cerr << "WARNING: corrupt line!" << endl;
            name.clear();
            content.clear();
            continue;
        }
        if( line[0] == '>' ){
            // Output complete entry
            if( !name.empty() ) {
                std::transform(content.begin(), content.end(), content.begin(), ::toupper);
                fastas.emplace_back( FASTA{name, content} );
            }
            // Save new entries header line
            if ( !line.empty() ) {
                name = line.substr(1);
            }
            // Clear content for new entry
            content.clear();
            // Found sequence for previous name
        } else if ( !name.empty() ){
            content += line;
        }
    }

    // Reached the end, so save final entry
    if ( !name.empty() ){
        std::transform(content.begin(), content.end(), content.begin(), ::toupper);
        fastas.emplace_back( name, content );
    }
}


void pfastq(const string& path, vector<FASTQ> &fastqs){
    /*
     * Parses FASTQ file and loads vector with FASTQ objects
     */
    ifstream file(path,
                  ios_base::in | ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    // Add gzip decompressor if needed
    if (path.substr(path.find_last_of('.') + 1) == "gz"){
        in.push(boost::iostreams::gzip_decompressor());
    }
    in.push(file);
    std::istream input{&in};
	if (!input.good()){
	    std::cerr << "Error opening '" << path << "'. Bailing out." << std::endl;
	}
	std::string line{}, name{}, seq{}, comment{}, quality{};
	while( std::getline( input, line ) ){
	    if ( line.empty() ){
	        cerr << "FASTQ contains empty line!" << endl;
	        continue;
	    }
	    if ( line[0] == '@' ){
	        name = line.substr(1);
	        std::getline( input, seq );
	        std::getline( input, line );
	        comment = line.substr(1);
	        std::getline(input, quality);

	        if( seq.size() == quality.size() ) {
                fastqs.emplace_back(name, seq, comment, quality);
            } else{
	            std::cout << "Read " << name << " is corrupt\n";
	        }
	    }
	}
}


void mmfastq(const string& path, vector<FASTQ>& fastqs){
    /*
     * Parses FASTQ file and loads vector with FASTQ objects
     */
    bio::stream<bio::mapped_file_source> file;
    file.open(bio::mapped_file_source(path));

    bio::filtering_streambuf<bio::input> in;

    // Add gzip decompressor if needed
    if (path.substr(path.find_last_of('.') + 1) == "gz"){
        in.push(bio::gzip_decompressor());
    }

    in.push(file);
    std::istream input{&in};
    if (!input.good()){
        std::cerr << "Error opening '" << path << "'. Bailing out." << std::endl;
    }

    std::string line{}, name{}, seq{}, comment{}, quality{};
    while( std::getline( input, line ) ){
        if ( line.empty() ){
            cerr << "FASTQ contains empty line!" << endl;
            continue;
        }
        if ( line[0] == '@' ){
            name = line.substr(1);
            std::getline( input, seq );
            std::getline( input, line );
            comment = line.substr(1);
            std::getline(input, quality);

            if( seq.size() == quality.size() ) {
                fastqs.emplace_back(name, seq, comment, quality);
            } else{
                std::cout << "Read " << name << " is corrupt\n";
            }
        }
    }
}
