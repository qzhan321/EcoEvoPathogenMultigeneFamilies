#ifndef varmodel_hpp
#define varmodel_hpp

#include "Population.hpp"

namespace varmodel {

void verify_simulation_state();


void run(bool override_seed, uint64_t random_seed);

}; // namespace varmodel

#endif // #ifndef varmodel_hpp
