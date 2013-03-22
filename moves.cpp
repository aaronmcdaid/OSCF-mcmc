#include "moves.hpp"

#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>
using namespace std;

#include "macros.hpp"

gsl_rng * r = NULL;
void			seed_the_random_number_generator(int seed) {
					assert( r == NULL );
					r = gsl_rng_alloc (gsl_rng_taus);
					gsl_rng_set(r, seed);
}


static vector<bool>	bernoullis_not_all_failed(const vector<long double> p_k) {
					const int64_t K = p_k.size();
					int num_of_successes = 0;
					int64_t attempts = 0;
					vector<bool> bools(K);
					while(num_of_successes == 0) {
						++attempts;
						for(int k = 0; k<K; ++k) {
							const bool b = gsl_ran_bernoulli(r, p_k.at(k));
							if(b)
								++ num_of_successes;
							bools.at(k) = b;
						}
					}
					// PP2(attempts, num_of_successes);
					return bools;
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
	{
		while( ! sc.state.get_edge_to_set_of_comms().at(e).empty() ) {
			int64_t comm_id_to_remove = * sc.state.get_edge_to_set_of_comms().at(e).begin();
			const OneCommunitySummary old_one_comm = sc.state.get_one_community_summary(comm_id_to_remove);
			sc.remove_edge(e, comm_id_to_remove);
			const OneCommunitySummary new_one_comm = sc.state.get_one_community_summary(comm_id_to_remove);
			delta_in_gibbs += sc.f(new_one_comm) - sc.f(old_one_comm);
		}
	}

	// For each of the K commmunities, how does the corresponding f change with the (re)addition of this edge?
	vector<long double> p_k(K);
	{
		for(int k = 0; k<K; ++k) {
			// const long double not_in = sc.score();
			const OneCommunitySummary old_one_comm = sc.state.get_one_community_summary(k);
			sc.add_edge(e, k);
			const OneCommunitySummary new_one_comm = sc.state.get_one_community_summary(k);
			const long double delta_score_one_edge = sc.f(new_one_comm) - sc.f(old_one_comm);
			// const long double is_in = sc.score();
			sc.remove_edge(e, k);
			// assert(not_in == sc.score());
			// assertVERYCLOSE(delta_score_one_edge, is_in - not_in);

			// Note: it *is* possible for is_in to be greater than not_in

			const long double extra_if_in = delta_score_one_edge;
			const long double a = exp2l(extra_if_in);
			const long double p = a / (1+a);
			// PP3(extra_if_in, a, p);
			assert(p > 0 && p < 1);
			p_k.at(k) = p;
			// PP2(k, p);
		}

	}

	// Assign the new values
	{
		const vector<bool> new_values_for_this_edge = bernoullis_not_all_failed(p_k);
		for(int k = 0; k<K; ++k) {
			if(new_values_for_this_edge.at(k)) {
				const OneCommunitySummary old_one_comm = sc.state.get_one_community_summary(k);
				sc.add_edge(e, k);
				const OneCommunitySummary new_one_comm = sc.state.get_one_community_summary(k);
				delta_in_gibbs += sc.f(new_one_comm) - sc.f(old_one_comm);
			}
		}
	}
	return delta_in_gibbs;
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
