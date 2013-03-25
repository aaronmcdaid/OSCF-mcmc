#ifndef SCORE_HPP__
#define SCORE_HPP__

#include "state.hpp"

#include<utility>
#include<tr1/unordered_map>
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

struct Score { // every modification *should* go through here eventually, so as to track the score.
	// But for now, this is just a passive object that recalculates all the scores from scratch each time
	State &		state;
	typedef std :: tr1 :: unordered_map < std :: pair<int64_t, int64_t>, std :: pair<long double, int64_t> > Cache_T;
	mutable Cache_T cache;
	explicit	Score(State & state_)	;
	long double	score()			const;
	long double	prior_on_K()		const;
	long double	product_on_fs()		const;
	long double	f(const int64_t num_edges, const int64_t num_unique_nodes_in_this_community)	const;
	long double	f(const OneCommunitySummary ocs)						const;
	long double 	what_would_change_if_I_deleted_an_empty_community()				const;

	// send the modifications through here

	long double	add_edge(int64_t e, int64_t comm_id);
	long double	remove_edge(int64_t e, int64_t comm_id);
	long double	append_empty_cluster()	;
	long double	delete_empty_cluster_from_the_end()	;
};

#endif
