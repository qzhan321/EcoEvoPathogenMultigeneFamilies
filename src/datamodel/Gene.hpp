#ifndef Gene_hpp
#define Gene_hpp

#include <array>
#include "parameters.hpp"

namespace varmodel {

enum GeneSource {
    SOURCE_POOL_ORIGINAL = 0,
    SOURCE_RECOMBINATION = 1,
    SOURCE_MUTATION = 2,
    SOURCE_POOL_MUTATION = 3,
};

struct Gene {
    Gene(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    uint64_t source;
    double is_functional;
    bool in_pool;
    uint64_t upsGroup;
    std::array<uint64_t, N_LOCI> alleles;
    double dur;
    double avTransRate;
    double originTime;
    double deathTime;
    double recombRate;
    uint64_t refcount;
};

} // namespace varmodel

#endif // #ifndef Gene_hpp
