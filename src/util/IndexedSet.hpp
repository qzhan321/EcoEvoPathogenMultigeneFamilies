#ifndef IndexedSet_hpp
#define IndexedSet_hpp

#include <vector>
#include <unordered_map>
#include <cassert>
#include <cstddef>
using std::size_t;

namespace varmodel {

template<typename T>
struct IndexedSet {
    std::vector<T *> vec;
    std::unordered_map<uint64_t, size_t> id_index_map; 
    
    void add(T * obj) {
        assert(id_index_map.find(obj->id) == id_index_map.end());
        vec.push_back(obj);
        id_index_map[obj->id] = vec.size() - 1;
    }
    
    void remove(T * obj) {
        size_t index = id_index_map[obj->id];
        remove_at_index(index);
    }
    
    T * remove_at_index(uint64_t index) {
        T * obj = vec[index];
        id_index_map.erase(obj->id);
        if(index < vec.size() - 1) {
            vec[index] = vec.back();
            id_index_map[vec[index]->id] = index;
        }
        vec.pop_back();
        return obj;
    }
    
    bool contains_id(uint64_t id) {
        return id_index_map.find(id) != id_index_map.end();
    }
    
    bool contains_object(T * obj) {
        return contains_id(obj->id);
    }
    
    T * object_for_id(int64_t id) {
        return object_at_index(id_index_map[id]);
    }
    
    T * object_at_index(size_t index) {
        return vec[index];
    }
    
    size_t size() {
        return vec.size();
    }
    
    std::vector<T *> & as_vector() {
        return vec;
    }
};

}; // namespace varmodel

#endif // #ifndef IndexedSet_hpp
