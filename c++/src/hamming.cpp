#include<string>
#include<vector>
#include<unordered_set>
#include<assert.h>

using namespace std;


int hamming_dist(const string& str1, const string& str2){
	if (str1.length() != str2.length()){
	       throw std::runtime_error{ "Strings must be same length!" };
	}

	int count = 0;
	for (size_t i = 0; i < str1.length(); i++)
	{
		if (str1[i] != str2[i])
			count++;
	}	
	return count;
}


int min_ham(const string& str, const vector<string>& background, const int& minimum){
	int ham;
	int min = str.length();
	for (const auto & i : background){
		ham = hamming_dist(str, i);
		if (ham < min){
			min = ham;
			// We know that there isn't a perfect match, so the lowest is actually 1
		} else if (min <= minimum){
		    break;
		}
	}
    assert( min < str.length() );
	return min;
}


int min_ham(const string& str,
            const vector<string>& background,
            const int& minimum,
            const int& start,
            const int& stop){
    int ham;
    int min = str.length();
    for (int i = start; i < stop; ++i){
        ham = hamming_dist(str, background[i]);
        if (ham < min){
            min = ham;
        } else if (min <= minimum){
            break;
        }
    }
    assert( min < str.length() );
    return min;
}
