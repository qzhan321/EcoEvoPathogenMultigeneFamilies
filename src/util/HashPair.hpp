#ifndef HashPair_hpp
#define HashPair_hpp

namespace varmodel {

template<typename T, typename Hash=std::hash<T>>
class HashPair
{
public:
	size_t operator()(std::pair<T,T> const & p) const
	{
		size_t hashFirst = _hash(p.first);
		hashFirst ^= _hash(p.second) + 0x9e3779b9 + (hashFirst << 6) + (hashFirst >> 2);
		return hashFirst;
	}
private:
	Hash _hash;
};

} // namespace varmodel

#endif // #ifndef HashPair_hpp
