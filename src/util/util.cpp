#include <cassert>
#include "util.hpp"
#include <sys/stat.h>

using namespace std;

namespace varmodel {

double add(std::vector<double> const & vec)
{
	double sum = 0.0;
	for(double val : vec) {
		sum += val;
	}
	return sum;
}

std::vector<double> addCumulative(std::vector<double> const & vec)
{
	vector<double> cumSum(vec);
	for(size_t i = 1; i < cumSum.size(); i++) {
		cumSum[i] += cumSum[i-1];
	}
	return cumSum;
}

bool file_exists(std::string const & filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

} // namespace varmodel
