#ifndef HASH_A_PAIR_HPP__
#define HASH_A_PAIR_HPP__
namespace std {
namespace tr1 {
	template<typename a, typename b>
	struct hash< std::pair<a, b> > {
		private:
		const hash<a> ah;
		const hash<b> bh;
		public:
		hash() : ah(), bh() {}
		size_t operator()(const std::pair<a, b> &p) const {
			return ah(-p.first) ^ bh(1+p.second);
		}
	};
}
}
#endif
