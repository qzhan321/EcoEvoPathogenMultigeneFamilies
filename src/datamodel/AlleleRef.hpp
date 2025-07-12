#ifndef AlleleRef_hpp
#define AlleleRef_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <unordered_map>

namespace varmodel {

struct AlleleRef {
    AlleleRef(uint64_t id) : id(id) { }
     
    uint64_t const id;
    
    uint64_t locus;
    uint64_t allele;
    double originTime;
};

} // namespace varmodel }

#endif // #define Immunity_hpp
