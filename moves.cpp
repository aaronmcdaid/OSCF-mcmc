#include "moves.hpp"

#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>
#include <cstdlib>
#include <algorithm>
using namespace std;

#include "macros.hpp"

gsl_rng * r = NULL;
void			seed_the_random_number_generator(int seed) {
					assert( r == NULL );
					r = gsl_rng_alloc (gsl_rng_taus);
					gsl_rng_set(r, seed);

					srand(seed); // I *think* this will seed std::random_shuffle
}


static pair< vector<bool>, long double>	bernoullis_not_all_failed(const vector<long double> p_k) {
					const int64_t K = p_k.size();
					long double log2_product_of_accepted_probabilities = 0.0L;

					long double p_rest_all_zeros = 0.0L;
					for(int k=0; k<K; ++k) {
						p_rest_all_zeros += log2l(1.0L - p_k.at(k));
					}

					int num_of_successes = 0;
					vector<bool> bools(K);

					for(int k = 0; k<K; ++k) {
						long double p = p_k.at(k);
						p_rest_all_zeros -= log2l(1.0L - p);
						{ // IF num_of_successes == 0, then we need to do something special
						  // in order to condition on there finally being at least one success
							if(num_of_successes==0) {
								p /= 1.0L - exp2l(p_rest_all_zeros + log2l(1.0L-p));
							}
						}
						const bool b = gsl_ran_bernoulli(r, p);
						if(b)
							++ num_of_successes;
						bools.at(k) = b;
						log2_product_of_accepted_probabilities += b ? log2l(p) : log2l(1.0L-p);
					}
					assertVERYCLOSE(p_rest_all_zeros , 0.0L);
					assert(num_of_successes > 0);
					assert(isfinite(log2_product_of_accepted_probabilities));
					return make_pair(bools, log2_product_of_accepted_probabilities);
}

long double		remove_edge_from_all_its_communities(int64_t e, Score &sc) {
		long double delta_in_here = 0.0L;
		while( ! sc.state.get_edge_to_set_of_comms().at(e).empty() ) {
			int64_t comm_id_to_remove = * sc.state.get_edge_to_set_of_comms().at(e).begin();
			delta_in_here += sc.remove_edge(e, comm_id_to_remove);
		}
		return delta_in_here;
}
long double		remove_edge_from_one_community_if_present(int64_t e, Score &sc, const int64_t comm_id_to_remove) {
		if( sc.state.get_edge_to_set_of_comms().at(e).count(comm_id_to_remove) == 1 )
			return sc.remove_edge(e, comm_id_to_remove);
		else
			return 0.0L;
}

long double		calculate_p_based_on_the_log_ratio(const long double delta_score_one_edge) {
			const long double a = exp2l(delta_score_one_edge);
			const long double p = a / (1+a);
			// PP3(extra_if_in, a, p);
			assert(isfinite(p));
			assert(p > 0 && p < 1);
			return p;
}

long double 		gibbsUpdate(int64_t e, Score & sc) {
	const int64_t K = sc.state.get_K();
	// Does NOT change K
	//
	// 1. remove the edge from *all* communities
	// 2. reassign to each in turn, calculating the delta-score in each case, giving the probability for that Bernoulli
	// 3. draw from the Bernoullis, but conditioning that it must be assigned to at least one community.

	// First, Remove the edge from *all* its communities
	long double delta_in_gibbs = 0.0L;
	delta_in_gibbs += remove_edge_from_all_its_communities(e, sc);

	// For each of the K commmunities, how does the corresponding f change with the (re)addition of this edge?
	vector<long double> p_k(K);
	{
		for(int k = 0; k<K; ++k) {
			// const long double not_in = sc.score();
			const long double delta_score_one_edge = sc.add_edge(e, k);
			// const long double is_in = sc.score();
			sc.state.remove_edge(e, k);
			// assert(not_in == sc.score());
			// assertVERYCLOSE(delta_score_one_edge, is_in - not_in);

			// Note: it *is* possible for is_in to be greater than not_in

			p_k.at(k) = calculate_p_based_on_the_log_ratio(delta_score_one_edge);
			// PP2(k, p);
		}

	}
	assert(sc.state.get_edge_to_set_of_comms().at(e).empty());

	// Assign the new values
	{
		const pair< vector<bool>,long double > new_values_for_this_edge = bernoullis_not_all_failed(p_k);
		for(int k = 0; k<K; ++k) {
			if(new_values_for_this_edge.first.at(k)) {
				delta_in_gibbs += sc.add_edge(e, k);
			}
		}
	}
	return delta_in_gibbs;
}
long double 		gibbsUpdateJustTwoComms(int64_t e, Score & sc, const int64_t main_cluster, const int64_t secondary_cluster) {
	// const int64_t K = sc.state.get_K();
	assert(main_cluster != secondary_cluster);
	// This is used with the split() and merge() moves.
	// It's based on gibbsUpdate,
	// *BUT* it works only on two communities
	//
	// 1. remove the edge from *both* communities (if not already removed)
	// 2. reassign to each in turn, calculating the delta-score in each case, giving the probability for that Bernoulli
	// 3. draw from the Bernoullis, but conditioning that it must be assigned to at least one community.

	// First, Remove the edge from *both* of these communities (only if it is in them, of course)
	long double delta_in_gibbs = 0.0L;
	delta_in_gibbs += remove_edge_from_one_community_if_present(e, sc, main_cluster);
	delta_in_gibbs += remove_edge_from_one_community_if_present(e, sc, secondary_cluster);

	// For each of the 2 commmunities, how does the corresponding f change with the (re)addition of this edge?
	vector<long double> p_k(2);
	{
		for(int justIterateOverTwo = 0; justIterateOverTwo < 2; ++ justIterateOverTwo) {
			const int k = justIterateOverTwo == 0 ? main_cluster : secondary_cluster;
			const long double delta_score_one_edge = sc.add_edge(e, k);
			sc.state.remove_edge(e, k); // Put things back the way they were
			p_k.at(justIterateOverTwo) = calculate_p_based_on_the_log_ratio(delta_score_one_edge);
		}

	}

	// Assign the new values
	{
		const pair< vector<bool>,long double > new_values_for_this_edge = bernoullis_not_all_failed(p_k);
		if(new_values_for_this_edge.first.at(0))
				delta_in_gibbs += sc.add_edge(e, main_cluster);
		if(new_values_for_this_edge.first.at(1))
				delta_in_gibbs += sc.add_edge(e, secondary_cluster);
	}
	return delta_in_gibbs;
}

struct TriState {
	short			x;
	explicit		TriState() : x(0)	{ }
	enum			Flags			{ IN_MAIN = 1, IN_SECN = 2 }; // Main, or Secondary, or both?
	void			put_in_MAIN()		{ x |= IN_MAIN; }
	void			put_in_SECN()		{ x |= IN_SECN; }
	bool			test_in_MAIN() const	{ return x & IN_MAIN; }
	bool			test_in_SECN() const	{ return x & IN_SECN; }
};

void			empty_both_clusters(
					const int main_cluster,
					const int secondary_cluster,
					const std :: tr1 :: unordered_map< int64_t , TriState > & original_state_of_these_edges,
					State & state
				) {
	For(edge_with_state, original_state_of_these_edges) {
		if(edge_with_state->second.test_in_MAIN()) { state.remove_edge(edge_with_state->first, main_cluster); }
		if(edge_with_state->second.test_in_SECN()) { state.remove_edge(edge_with_state->first, secondary_cluster); }
	}
	assert(state.get_comms().at(main_cluster     ).empty());
	assert(state.get_comms().at(secondary_cluster).empty());
}

long double		metroK(Score & sc) {
					// Either add or remove a cluster at random
					// This will affect the prior on K obviously, but don't forget the f(0,0) term
					if(gsl_ran_bernoulli(r, 0.5)) {
						// cout << "Attempt append" << endl;
						// Attempt to Add an empty cluster
						/// const long double pre = sc.score();
						const long double delta_score = sc.append_empty_cluster();
						/// const long double post = sc.score();
						assert(delta_score < 0);
						/// assert(VERYCLOSE(post - pre, delta_score));
						if( log2l(gsl_rng_uniform(r)) < delta_score ) {
							// Accept
							// .. but let's move it to a random location
							const int64_t target_cluster_id = sc.state.get_K() * gsl_rng_uniform(r);
							sc.state.swap_cluster_to_the_end(target_cluster_id);
							return delta_score;
						} else {
							// Reject
							sc.delete_empty_cluster_from_the_end();
							return 0.0L;
						}
					} else {
						// cout << "Attempt removal" << endl;
						// Attempt to Remove an empty cluster

						// First, select a cluster at random to be our target.
						// If it's empty, just bail out immediately
						assert(sc.state.get_K() > 0);
						const int64_t target_cluster_id = sc.state.get_K() * gsl_rng_uniform(r);
						if(!sc.state.get_comms().at(target_cluster_id).empty()) {
							return 0.0L;
						}

						// OK, it's empty , we can propose deleting it

						/// const long double pre_ = sc.score();
						sc.state.swap_cluster_to_the_end(target_cluster_id);
						/// const long double post_ = sc.score();
						/// assert(VERYCLOSE(pre_,post_));

						/// const long double pre = sc.score();
						const long double delta_score = sc.delete_empty_cluster_from_the_end();
						/// const long double post = sc.score();
						assert(delta_score > 0);
						/// assert(VERYCLOSE(post - pre, delta_score));
						if( log2l(gsl_rng_uniform(r)) < delta_score ) {
							// Accept
							return delta_score;
						} else {
							// Reject
							assert(1==2); // It should always Accept this proposal
							sc.append_empty_cluster();
							sc.state.swap_cluster_to_the_end(target_cluster_id);
							return 0.0L;
						}
					}
}
