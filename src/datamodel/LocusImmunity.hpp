#ifndef LocusImmunity_hpp
#define LocusImmunity_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <unordered_map>

namespace varmodel {

struct LocusImmunity {
    LocusImmunity(uint64_t id) : id(id) { }
     
    uint64_t const id;
    std::unordered_map<uint64_t,double> immunity_level_by_allele;
    
};

} // namespace varmodel }

#endif // #define Immunity_hpp
