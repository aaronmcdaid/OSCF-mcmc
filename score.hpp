#include "state.hpp"
struct Score { // every modification *should* go through here eventually, so as to track the score.
	// But for now, this is just a passive object that recalculates all the scores from scratch each time
	State &		state;
	explicit	Score(State & state_)	: state(state_) {
	}
	long double	score()			const;
	long double	prior_on_K()		const;
	long double	product_on_fs()		const;
	long double	f(const int64_t num_edges, const int64_t num_unique_nodes_in_this_community)	const;
};

