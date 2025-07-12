#ifndef random_hpp
#define random_hpp

#include <vector>
#include <unordered_set>
#include <random>
#include <cassert>
#include "util.hpp"
#include "HashPair.hpp"

namespace varmodel {

typedef std::mt19937_64 rng_t;

template<typename T>
T drawUniformIndex(rng_t & rng, T size) {
	std::uniform_int_distribution<T> unif(0, size - 1);
	return unif(rng);
}

template<typename T>
T drawUniformIndexExcept(rng_t & rng, T size, T except) {
	std::uniform_int_distribution<T> unif(0, size - 2);
	T value = unif(rng);
	if(value >= except) {
		value++;
	}
	return value;
}

template<typename T>
std::vector<T> drawUniformIndicesExcept(rng_t & rng, T size, T count, T except) {
	std::uniform_int_distribution<T> unif(0, size - 2);
	std::vector<T> indexes(count);
	for(T i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			if(indexes[i] >= except) {
				indexes[i]++;
			}
			done = true;
			for(T j = 0; j < i; j++) {
				if(indexes[i] == indexes[j]) {
					done = false;
					break;
				}
			}
		}
	}
	sort(indexes.begin(), indexes.end());
	return indexes;
}

template<typename T>
T drawUniformIndexExcept(rng_t & rng, T size, std::vector<T> except) {
	sort(except.begin(), except.end());
	
	std::uniform_int_distribution<T> unif(0, size - 1 - T(except.size()));
	T value = unif(rng);
	for(T exceptVal : except) {
		if(value >= exceptVal) {
			value++;
		}
	}
	return value;
}

template<typename T>
std::vector<T> drawUniformIndicesExcept(rng_t & rng, T size, T count, std::vector<T> except) {
	sort(except.begin(), except.end());
	
	std::uniform_int_distribution<T> unif(0, size - 2);
	std::vector<T> indexes(count);
	for(T i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			for(T exceptVal : except) {
				if(indexes[i] >= exceptVal) {
					indexes[i]++;
				}
			}
			done = true;
			for(T j = 0; j < i; j++) {
				if(indexes[i] == indexes[j]) {
					done = false;
					break;
				}
			}
		}
	}
	sort(indexes.begin(), indexes.end());
	return indexes;
}

template<typename T>
std::vector<T> drawUniformIndices(rng_t & rng, T size, T count, bool sorted) {
	std::uniform_int_distribution<T> unif(0, size - 1);
	std::vector<T> indexes(count);
	std::unordered_set<T> indexSet;
	for(T i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			if(indexSet.find(indexes[i]) != indexSet.end()) {
				done = false;
			}
			else {
				indexSet.insert(indexes[i]);
				done = true;
			}
		}
	}
	if(sorted) {
		sort(indexes.begin(), indexes.end());
	}
	return indexes;
}

template<typename T>
std::vector<T> drawRandomSubset(rng_t & rng, std::vector<T> srcVec, size_t count, bool sortedByIndex)
{
	std::vector<T> dstVec;
	dstVec.reserve(count);
	for(size_t index : drawUniformIndices(rng, srcVec.size(), count, sortedByIndex)) {
		dstVec.push_back(srcVec[index]);
	}
	return dstVec;
}

template<typename T>
std::vector<T> drawUniformIndices(rng_t & rng, T size, T count) {
	return drawUniformIndices(rng, size, count, true);
}

template<typename T>
std::vector<T> drawMultipleBernoulli(rng_t & rng, T size, double p, bool sortByIndex) {
	std::binomial_distribution<int64_t> countDist(size, p);
	int64_t count = countDist(rng);
	std::vector<T> indices = drawUniformIndices(rng, size, T(count), sortByIndex);
	assert(indices.size() == count);
	return indices;
}

template<typename T>
std::vector<T> drawMultipleBernoulli(rng_t & rng, T size, double p) {
	return drawMultipleBernoulli(rng, size, p, true);
}

template<typename T>
std::vector<T> drawMultipleBernoulliRandomSubset(rng_t & rng, std::vector<T> srcVec, double p, bool sortedByIndex)
{
	std::vector<T> dstVec;
	std::vector<size_t> indices = drawMultipleBernoulli(rng, srcVec.size(), p, sortedByIndex);
	dstVec.reserve(indices.size());
	for(size_t index : indices) {
		dstVec.push_back(srcVec[index]);
	}
	return dstVec;
}

template<typename T>
std::vector<std::pair<T, T>> drawUniformIndexPairs(rng_t & rng, T size, T count, bool ordered, bool different) {
	if(ordered && different) {
		assert(count <= size * (size - 1));
	}
	else if(ordered && !different) {
		assert(count <= size * size);
	}
	else if(!ordered && different) {
		assert(count <= size * (size - 1) / 2);
	}
	else {
		assert(count <= size * size);
	}
	
	std::uniform_int_distribution<T> unif(0, size - 1);
	std::vector<std::pair<T, T>> indexPairs(count);
	std::unordered_set<std::pair<T, T>, HashPair<T>> indexPairSet;
	for(T i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			T i1 = unif(rng);
			T i2 = unif(rng);
			if(!different || (i1 != i2)) {
				if(!ordered && (i2 < i1)) {
					std::swap(i1, i2);
				}
				indexPairs[i] = std::pair<T,T>(i1, i2);
				if(indexPairSet.find(indexPairs[i]) != indexPairSet.end()) {
					done = false;
				}
				else {
					indexPairSet.insert(indexPairs[i]);
					done = true;
				}
			}
		}
	}
	return indexPairs;
}

template<typename T>
std::vector<std::pair<T, T>> drawMultipleBernoulliIndexPairs(rng_t & rng, T size, double p, bool ordered, bool different) {
	T sizePairs;
	if(ordered && different) {
		sizePairs = size * (size - 1);
	}
	else if(ordered && !different) {
		sizePairs = size * size;
	}
	else if(!ordered && different) {
		sizePairs = size * (size - 1) / 2;
	}
	else {
		sizePairs = size * size;
	}
	
	std::binomial_distribution<int64_t> countDist(sizePairs, p);
	int64_t count = countDist(rng);
	std::vector<std::pair<T,T>> indices = drawUniformIndexPairs(rng, size, T(count), ordered, different);
	assert(indices.size() == count);
	return indices;
}

size_t sampleDiscreteLinearSearch(rng_t & rng, std::vector<double> const & weights);

} // namespace varmodel

#endif // #ifndef random_hpp
