#ifndef EventQueue_hpp
#define EventQueue_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <cassert>
#include <iostream>

namespace varmodel {

template<typename T, double (*V)(T *)>
bool Compare(T * o1, T * o2) {
    if(V(o1) == V(o2)) {
        return o1->id < o2->id;
    }
    return V(o1) < V(o2);
}

template<typename T, double (*V)(T *)>
class EventQueue
{
public:
	EventQueue():
		heap(0),
		indexes(0)
	{
	}
	
	EventQueue(std::vector<T> const & initVec):
		heap(initVec),
		indexes(makeIndexes())
	{
		buildHeap();
	}
	
	void update(T * obj)
	{
		assert(indexes.find(obj) != indexes.end());
		updateAtHeapIndex(indexes[obj]);
	}
	
	bool add(T * value)
	{
		if(indexes.find(value) != indexes.end()) {
			return false;
		}
		
		heap.push_back(value);
		size_t index = heap.size() - 1;
		indexes[value] = index;
		
		assert(indexes.size() == heap.size());
		heapifyUp(index);
		
		return true;
	}
	
	bool remove(T * value)
	{
		auto itr = indexes.find(value);
		if(itr == indexes.end()) {
			return false;
		}
		
		size_t index = itr->second;
		if(index == heap.size() - 1) {
			heap.pop_back();
			indexes.erase(itr);
		}
		else {
			T * obj = heap.back();
			heap.pop_back();
			heap[index] = obj;
			indexes.erase(itr);
			indexes[obj] = index;
			updateAtHeapIndex(index);
		}
		assert(indexes.size() == heap.size());
		
		return true;
	}
	
    double next_time() {
        if(heap.size() == 0) {
            return std::numeric_limits<double>::infinity();
        }
        return V(head());
    }
    
	T * head()
	{
		if(heap.size() == 0) {
			return NULL;
		}
		
		return heap[0];
	}
	
	size_t getMaxParentIndex()
	{
		assert(heap.size() > 1);
		return (heap.size() - 2)/2;
	}
	
	bool verify_heap()
	{
		if(heap.size() <= 1) {
			return true;
		}
		
		size_t maxParentIndex = getMaxParentIndex();
		for(size_t i = 0; i <= maxParentIndex; i++)
		{
			size_t left = 2*i + 1;
			assert(left < heap.size());
			if(Compare<T, V>(heap[left], heap[i])) {
				return false;
			}
			
			size_t right = left + 1;
			if(right < heap.size()) {
				if(Compare<T, V>(heap[right], heap[i])) {
					return false;
				}
			}
		}
		assert(2 * (maxParentIndex + 1) + 1 >= heap.size());
		return true;
	}
	
	size_t size()
	{
		return heap.size();
	}
	
	bool contains(T const & value)
	{
		return indexes.find(value) != indexes.end();
	}
	
private:
	std::vector<T *> heap;
	std::unordered_map<T *, size_t> indexes;
	
	EventQueue(std::vector<T *> const & initVec, bool useAsInitialHeap):
		heap(initVec),
		indexes(makeIndexes())
	{
		if(!useAsInitialHeap) {
			buildHeap();
		}
	}
	
	std::unordered_map<T *, size_t> makeIndexes()
	{
		std::unordered_map<T *, size_t> tmpIndexes(heap.size());
		for(size_t i = 0; i < heap.size(); i++)
		{
			assert(tmpIndexes.find(heap[i]) == tmpIndexes.end());
			indexes[heap[i]] = i;
		}
		
		return tmpIndexes;
	}
	
	void buildHeap()
	{
		if(heap.size() < 2) {
			return;
		}
		
		// Build heap
		size_t maxParentIndex = getMaxParentIndex();
		for(size_t i = maxParentIndex + 1; i > 0; i--) {
			heapifyDown(i - 1);
		}
	}
	
	bool heapifyUp(size_t i)
	{
		assert(i  < heap.size());
		
		bool moved = false;
		while(i > 0) {
			size_t parent = (i-1)/2;
			if(Compare<T, V>(heap[i], heap[parent])) {
				swap(i, parent);
				i = parent;
				moved = true;
			}
			else break;
		}
		return moved;
	}
	
	bool heapifyDown(size_t i)
	{
		return heapifyDown(i, false);
	}
	
	bool heapifyDown(size_t i, bool removing)
	{
		assert(i < heap.size());
		if(heap.size() < 2) {
			return false;
		}
		
		bool moved = false;
		
		size_t size = heap.size();
		size_t maxParentIndex = (size - 2) / 2;
		while(i <= maxParentIndex)
		{
			T * nVal = heap[i];
			
			size_t left = 2*i + 1;
			if(left >= size) break;
			T * leftVal = heap[left];
			
			size_t right = left + 1;
			T * rightVal = right < size ? heap[right] : NULL;
			
			size_t min = i;
			
			if(removing || Compare<T, V>(leftVal, nVal))
			{
				if(rightVal != NULL && Compare<T, V>(rightVal, leftVal)) {
					min = right;
				}
				else {
					min = left;
				}
			}
			else if(rightVal != NULL && Compare<T, V>(rightVal, nVal)) {
				min = right;
			}
			else {
				break;
			}
			
			moved = true;
			swap(i, min);
			i = min;
		}
		return moved;
	}
	
	void updateAtHeapIndex(size_t index)
	{
		if(!heapifyUp(index))
		{
			heapifyDown(index);
		}
	}
	
	void swap(size_t i1, size_t i2)
	{
		assert(i1 < heap.size());
		assert(i2 < heap.size());
		assert(i1 != i2);
		
		T * e1 = heap[i2];
		T * e2 = heap[i1];
		heap[i1] = e1;
		heap[i2] = e2;
		
		indexes[e1] = i1;
		indexes[e2] = i2;
	}
};

} // namespace varmodel

#endif // #ifndef EventQueue_hpp
