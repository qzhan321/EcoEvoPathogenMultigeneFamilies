#ifndef Host_hpp
#define Host_hpp

#include <stdint.h>
#include <sqlite3.h>
#include <unordered_set>

namespace varmodel {

struct Population;
struct Infection;
struct ImmuneHistory;

struct Host {
    Host(uint64_t id) : id(id) { }
    
    uint64_t const id;
    
    Population * population;
    
    double birth_time;
    double death_time;
    
    double next_immunity_loss_time;
    uint64_t completed_infection_count;
    
    // One-to-many relationships
    std::unordered_set<Infection *> infections;
    
    ImmuneHistory * immune_history;
    
    bool MDA_effective_period;
};

} // namespace varmodel

#endif // #ifndef Host_hpp
