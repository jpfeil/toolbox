// Library for kmer counting tasks
#include<algorithm>
#include<cassert>
#include<fstream>
#include<iostream>
#include<stdexcept>
#include<string>
#include<vector>
#include<unordered_map>
#include<cstdio>
#include<unordered_set>
#include<string_view>

#include<boost/iostreams/device/mapped_file.hpp>
#include<boost/iostreams/filtering_streambuf.hpp>
#include<boost/iostreams/filter/gzip.hpp>
#include<boost/interprocess/file_mapping.hpp>
#include<boost/iostreams/stream.hpp>


#include<vaccinaTE/parser.h>
#include<vaccinaTE/util.h>
#include<vaccinaTE/kmer.h>

namespace bip = boost::interprocess;
namespace bio = boost::iostreams;


int count_lines(const std::string& path){
    char command[150];
    char buffer[150];
    std::string slines;
    if (path.substr(path.find_last_of('.') + 1) == "gz"){
        sprintf(command, "zcat %s | wc -l", path.c_str());
    } else{
        sprintf(command, "cat %s | wc -l", path.c_str());
    }

    FILE* pipe = popen( command, "r" );
    if (pipe){
        while ( !feof( pipe )){
            if ( fgets( buffer, 150, pipe ) != nullptr ){
                slines.append( buffer );
            }
        }
        pclose(pipe);
    }

    int lines = 0;
    try{
        lines = std::stoi(slines);

    } catch( std::invalid_argument ){
        std::cout << "Could not infer lines: " << slines << '\n';
    }

    return lines;
}


void populate(const int& length,
              const std::string& sequence,
              unordered_set<std::string> &container){
    /*
     * Takes a string and returns a vector full of kmers
     */
	for (size_t i = 0; i < sequence.size() - length + 1; i++){
	    if ( i + length <= sequence.size() ){
            container.insert( sequence.substr(i, length) );
	    }
	}
}


unordered_set<std::string> batch_populate(const int& length, const vector<FASTA>& sequences){
    /*
     * Takes a string and returns a vector full of kmers
     */
    unordered_set<std::string> container{};
    for (const auto& sequence: sequences){
        populate(length, sequence.seq, container);
    }
    return container;
}


void kmer_map(const std::string& name,
              const int& length,
              const std::string& sequence,
              unordered_map<std::string, vector<std::string>>& kmap){
    /*
     * Takes a string and places all kmers into map
     * data structure.
     */
    std::string kmer{};
    for (size_t i=0; i < sequence.size() - length + 1; i++){
        kmer = sequence.substr(i, length);
        if(kmap.find(kmer) == kmap.end() ){
            kmap[kmer] = {name};
        } else if (kmap.find(kmer) != kmap.end() ){
            if (find(kmap[kmer].begin(), kmap[kmer].end(), name) == kmap[kmer].end()){
                kmap[kmer].push_back(name);
            }
        } else {
            throw invalid_argument( "Kmer map is corrupt!" );
        }
    }
}


void kmer_count(const int& length,
                const std::string& sequence,
                std::unordered_map<std::string, int>& counter){
	std::string kmer{};
	for (size_t i=0; i < sequence.size() - length + 1; i++){
		kmer = sequence.substr(i, length);
		if( counter.find(kmer) ==  counter.end() ){
			counter[kmer] = 1;	
		} else if ( counter.find(kmer) != counter.end()){
			counter[kmer] = counter[kmer] + 1;
		} else {
			throw invalid_argument( "Kmer map is corrupt!" );
		}
	}
}


void kmer_finder(const std::string& read,
                 const std::string& quality,
                 const std::string& readtype,
                 const std::string& libtype,
                 const int& min_quality,
                 Counter& counter) {

    // Create a vector in case it's a IU
    // library and we need to check both
    // strands for a match.
    std::vector<std::string> reads{};
    std::vector<std::string> qualities{};

    // Copy the quality and reverse it.
    std::string rev_quality = quality;
    reverse(rev_quality.begin(), rev_quality.end());

    // Read 1 is reverse complemented to the
    // sense strand for the SR and ISR library types
    if (readtype == "R1" && (libtype == "SR" | libtype == "ISR")) {
        reads.emplace_back( revcomp(read) );
        qualities.emplace_back( rev_quality );

        // Read 2 is reverse complemented to the
        // sense strand for the SF and ISF library types
    } else if (readtype == "R2" && (libtype == "SF" | libtype == "ISF")) {
            reads.emplace_back( revcomp(read) );
            qualities.emplace_back(rev_quality);

    } else if ( libtype == "IU "){
        reads.emplace_back( read );
        qualities.emplace_back( quality );

        reads.emplace_back( revcomp(read) );
        qualities.emplace_back( rev_quality );

        // Otherwise just save the string
    } else {
        reads.emplace_back( read );
        qualities.emplace_back( quality );
    }

    // Generate all kmer lengths and increment
    // a counter if the kmer matches one in the counter
    std::string kmer{};
    std::string_view kmer_qual{};
    bool keep;

    for (size_t r = 0; r < reads.size(); r++ ){
        std::string_view qview{qualities[r]};
        for (size_t i = 0; i < reads[r].size(); i++) {
            for (size_t k: counter.kmer_lengths) {
                keep = true;

                // If the next kmer is longer than the
                // read continue.
                if (i + k > reads[r].size()){
                    break;
                }

                // Get the read kmer and the quality substring
                kmer = reads[r].substr(i, k);
                kmer_qual = qview.substr(i, k);

                // Check for minimum base calling score
                for (auto& score: kmer_qual){
                    if ( int(score) < min_quality ){
                        keep = false;
                        break;
                    }
                }

                if (keep) {
                    counter.increment(kmer);
                }
            }
        }
    }
}


void count_reads(const std::string& path,
                 const std::string& libtype,
                 const std::string& readtype,
                 Counter& counter,
                 int& total,
                 const int& min_quality) {
    bio::stream<bio::mapped_file_source> file;
    file.open(bio::mapped_file_source(path));

    bio::filtering_streambuf<bio::input> in;

    std::string_view path_view{path.c_str(), path.size()};

    // Add gzip decompressor if needed
    if (path_view.substr(path_view.find_last_of('.') + 1) == "gz") {
        in.push(bio::gzip_decompressor());
    }

    in.push(file);
    std::istream input{&in};
    if (!input.good()) {
        std::cerr << "Error opening '" << path << "'. Bailing out." << '\n';
    }

    std::string line{}, read{}, quality{};

    while (std::getline(input, line)) {
        // Read in header
        assert(line[0] == '@');

        // Read sequence line
        std::getline(input, read);
        total += 1;

        // Read comment line
        std::getline(input, line);

        // Read quality line
        std::getline(input, quality);
        assert(read.size() == quality.size());

        kmer_finder(read, quality, readtype, libtype, min_quality, counter);
    }
}


void probe_reads(const std::string &path,
                 const std::string &libtype,
                 const std::string &readtype,
                 Counter &vcounter,
                 Counter &mcounter,
                 Counter &ncounter,
                 int &total,
                 const int &min_quality) {
    bio::stream<bio::mapped_file_source> file;
    file.open(bio::mapped_file_source(path));
    bio::filtering_streambuf<bio::input> in;

    std::string_view path_view{path.c_str(), path.size()};

    // Add gzip decompressor if needed
    if (path_view.substr(path_view.find_last_of('.') + 1) == "gz") {
        in.push(bio::gzip_decompressor());
    }
    in.push(file);
    std::istream input{&in};
    if (!input.good()) {
        std::cerr << "Error opening '" << path << "'. Bailing out." << '\n';
    }
    std::string line{}, name{}, read{}, quality{};
    while (std::getline(input, line)) {
        assert(line[0] == '@');
        std::string_view header_view{line.c_str(), line.size()};
        name = header_view.substr(1);

        // Read sequence line
        std::getline(input, read);

        // Read comment line
        std::getline(input, line);

        // Read quality line
        std::getline(input, quality);
        assert(read.size() == quality.size());

        kmer_finder(read, quality, readtype, libtype, min_quality, vcounter);
        kmer_finder(read, quality, readtype, libtype, min_quality, mcounter);
        kmer_finder(read, quality, readtype, libtype, min_quality, ncounter);
    }
}

