#ifndef ImmuneHistory_hpp
#define ImmuneHistory_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <array>
#include "parameters.hpp"
#include "IndexedSet.hpp"

namespace varmodel {

struct LocusImmunity;
struct AlleleRef;

struct ImmuneHistory {
    ImmuneHistory(uint64_t id) : id(id) { }
     
    uint64_t const id;
    
    // Immunity levels at each locus
    std::array<LocusImmunity *, N_LOCI> immunity_by_locus;
    
    // Random-samplable set of references to each allele immunity currently present
    IndexedSet<AlleleRef> alleles_with_immunity; 
};

} // namespace varmodel }

#endif // #define ImmuneHistory_hpp
