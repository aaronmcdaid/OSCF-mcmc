#ifndef SCORE_HPP__
#define SCORE_HPP__

#include "state.hpp"

#include<vector>
#include<utility>
#include<tr1/unordered_map>
#include"hash_a_pair.hpp"

struct Score { // every modification *should* go through here eventually, so as to track the score.
	// But for now, this is just a passive object that recalculates all the scores from scratch each time
	State &		state;
	typedef std :: tr1 :: unordered_map < int64_t , std :: pair<long double, int64_t> > Cache_T;
	mutable std :: vector<Cache_T> cache; // One  cache for each number between 0 and N (inclusive). i.e. N+1 entries
	explicit	Score(State & state_)	;
	long double	score()			const;
	long double	prior_on_K()		const;
	long double	product_on_fs()		const;
	long double	f(const int64_t num_edges, const int64_t num_unique_nodes_in_this_community)	const;
	long double	f(const OneCommunitySummary ocs)						const;
	long double 	what_would_change_if_I_deleted_an_empty_community()				const;
	long double 	what_would_change_if_I_added_an_empty_community()				const;

	// send the modifications through here

	long double	add_edge(int64_t e, int64_t comm_id);
	long double	remove_edge(int64_t e, int64_t comm_id);
	long double	add_edge_if_not_already(int64_t e, int64_t comm_id);
	long double	remove_edge_if_not_already(int64_t e, int64_t comm_id);
	long double	set(int64_t e, int64_t comm_id, bool b);
	long double	append_empty_cluster()	;
	long double	delete_empty_cluster_from_the_end()	;

	// include some query funcs here, so they can be inlined
	long double	if_this_edge_is_added(const int64_t comm, const int64_t edge) const {
		const network :: EdgeSet :: Edge edge_details = this->state.net->edge_set->edges.at(edge);
		return this->if_this_edge_is_added(comm, edge_details.left, edge_details.right);
	}
	long double	if_this_edge_is_added(const int64_t comm, const int64_t left, const int64_t right) const {
			const Community & that_comm = this->state.get_comms().at(comm);
			const bool left_in_that_comm_already = that_comm.test_node(left);
			const bool right_in_that_comm_already = that_comm.test_node(right);
			assert(left != right);
			const int num_extra_nodes = (left_in_that_comm_already?0:1) + (right_in_that_comm_already?0:1);

			const int num_nodes_pre = that_comm.get_num_unique_nodes_in_this_community();
			const int num_nodes_post = that_comm.get_num_unique_nodes_in_this_community() + num_extra_nodes;

			const int64_t num_edges_pre = that_comm.get_num_edges();
			const int64_t num_edges_post = num_edges_pre + 1;
			const long double expected_delta_score = this->f(num_edges_post, num_nodes_post) - this->f(num_edges_pre, num_nodes_pre);
			return expected_delta_score;
	}
};

#endif
