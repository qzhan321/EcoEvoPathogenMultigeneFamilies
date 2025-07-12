#include "random.hpp"

using namespace std;

#include <vector>
#include <iostream>
#include <cassert>
#include "util.hpp"

namespace varmodel {

size_t sampleDiscreteLinearSearch(rng_t & rng, std::vector<double> const & weights)
{
	size_t n = weights.size();
	
	assert(n > 0);
	if(n == 1) {
		return 0;
	}
	
	vector<double> cumSum = addCumulative(weights);
	uniform_real_distribution<double> ud(0, cumSum[n - 1]);
	double u = ud(rng);
	for(size_t i = 0; i < n; i++) {
		if(u < cumSum[i]) {
			return i;
		}
	}
	return n - 1;
}

} // namespace varmodel
