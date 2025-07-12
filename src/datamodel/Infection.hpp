#ifndef Infection_hpp
#define Infection_hpp

#include <stdint.h>
#include <sqlite3.h>

#include "parameters.hpp"
#include <array>

namespace varmodel {

struct Strain;
struct Host;
struct Gene;

struct Infection { 
    Infection(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Strain * strain;
    Host * host;
    
    int64_t expression_index;
    
    double transition_time;
    double mutation_time;
    double recombination_time;
    double clearance_time;
    double infected_time;
    double totalRecombRate;

    std::array<Gene *, N_GENES_PER_STRAIN> expression_order;
    std::array<double, N_GENES_PER_STRAIN*(N_GENES_PER_STRAIN-1)/2> pair_recomb_rates;
};

} // namespace varmodel

#endif // #ifndef Infection_hpp
