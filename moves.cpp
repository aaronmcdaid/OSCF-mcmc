#include "moves.hpp"

#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <algorithm>
using namespace std;

#include "macros.hpp"

gsl_rng * r = NULL;

static long double probability_of_selecting_these_two_comms(const int main_cluster, const int secondary_cluster, const State &st);

void			seed_the_random_number_generator(int seed) {
					assert( r == NULL );
					r = gsl_rng_alloc (gsl_rng_taus);
					gsl_rng_set(r, seed);

					srand(seed); // I *think* this will seed std::random_shuffle
}


static pair< vector<bool>, long double>	bernoullis_not_all_failed(
							const vector<long double> p_k
							, std :: pair<int,int> possibly_force = std :: make_pair(-1,-1)
							) {
					const int64_t K = p_k.size();
	{
		if(possibly_force.first != -1) {
			assert(possibly_force.first  >= 0 && possibly_force.first  < 2); // It's just zero or one
			assert(possibly_force.second >= 0 && possibly_force.second < 2); // It's just zero or one
			assert(K==2);
		}
	}
					long double log2_product_of_accepted_probabilities = 0.0L;
#define log2_one_plus_l(x) (M_LOG2E * log1pl(x))

					long double _p_rest_all_zeros = 0.0L;
					vector<long double> p_rest_all_zeros_vector(K+1);
					p_rest_all_zeros_vector.at(K) = 0.0L;
					for(int k=K-1; k>=0; --k) {
						         p_rest_all_zeros_vector.at(k) = p_rest_all_zeros_vector.at(k+1) + log2_one_plus_l(- p_k.at(k));
						assertVERYCLOSE(p_rest_all_zeros_vector.at(k) , p_rest_all_zeros_vector.at(k+1) + log2l(1.0L - p_k.at(k))       );
						assert(isfinite(p_rest_all_zeros_vector.at(k)));
						_p_rest_all_zeros += log2l(1.0L - p_k.at(k));
						//PP4(k, p_k.at(k), _p_rest_all_zeros, p_rest_all_zeros_vector.at(k));
					}
					// for(int k=0; k<K; ++k) { PP2(k, p_rest_all_zeros_vector.at(k)); }
					assert(isfinite(_p_rest_all_zeros));

					int num_of_successes = 0;
					vector<bool> bools(K);
					for(int k = 0; k<K; ++k) {
						// - Take the unconditional p
						// - Calculate the unconditional ratio
						// - Calculate the   conditional ratio
						// - It should still be a finate number

						const long double NEW_p_rest_all_zeros = p_rest_all_zeros_vector.at(k+1);

						const long double uncond_p = p_k.at(k);
						assert(uncond_p > 0 && uncond_p < 1);

						const long double log2_uncond_ratio = log2l(uncond_p) - log2_one_plus_l(-uncond_p);
						assert(isfinite(log2_uncond_ratio));

						long double log2_cond_ratio = log2_uncond_ratio;
						if(num_of_successes == 0) {
							long double shrink_for_the_conditionality = - log2_one_plus_l(-exp2l(NEW_p_rest_all_zeros));
							if(k+1==K)
								assertVERYCLOSE(NEW_p_rest_all_zeros, 0.0L);
							if(isinfl(shrink_for_the_conditionality) == 1) {
								shrink_for_the_conditionality = std :: numeric_limits<long double> :: max(); // to effectively guarantee a success
							}
							assert(shrink_for_the_conditionality > 0);
							assert(isfinite(shrink_for_the_conditionality));
							log2_cond_ratio += shrink_for_the_conditionality;
						}
						assert(isfinite(log2_cond_ratio));
						assert(log2_cond_ratio >= log2_uncond_ratio);

						long double cond_p_via_ratio = exp2l( - log2_one_plus_l( exp2l(-log2_cond_ratio) ) ) ;
						if(k+1==K && num_of_successes==0)
							cond_p_via_ratio = 1.0L; // must be forced on if it's the last one
						assert(isfinite(cond_p_via_ratio));
						assert(cond_p_via_ratio > 0);
						assert(cond_p_via_ratio <= 1);
						bool b = false;
						if(possibly_force.first!=-1) {
							if(k==0) b = possibly_force.first;
							if(k==1) b = possibly_force.second;
						} else
							b = gsl_ran_bernoulli(r, cond_p_via_ratio);
						if(b)
							++ num_of_successes;
						bools.at(k) = b;
						assert(isfinite(log2_product_of_accepted_probabilities));
						if(b==0 && cond_p_via_ratio == 1.0L) { // Can only happen with possibly_force
							log2_product_of_accepted_probabilities += - std :: numeric_limits<long double> :: max();
						} else {
							log2_product_of_accepted_probabilities += b ? log2l(cond_p_via_ratio) : log2_one_plus_l(-cond_p_via_ratio);
						}
						assert(isfinite(log2_product_of_accepted_probabilities));
					}
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
pair<long double,long double> 	gibbsUpdateJustTwoComms(
							int64_t e
							, Score & sc, const int64_t main_cluster
							, const int64_t secondary_cluster
							, pair<int,int> possibly_force /* DEFAULTED = make_pair(-1,-1)*/ 
							)
{
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
	const pair< vector<bool>,long double > new_values_for_this_edge = bernoullis_not_all_failed(p_k, possibly_force);
	{
		if(new_values_for_this_edge.first.at(0))
				delta_in_gibbs += sc.add_edge(e, main_cluster);
		if(new_values_for_this_edge.first.at(1))
				delta_in_gibbs += sc.add_edge(e, secondary_cluster);
	}
	return make_pair(delta_in_gibbs, new_values_for_this_edge.second); // includes the proposal probability
}

struct TriState {
	short			x;
	explicit		TriState() : x(0)	{ }
	enum			Flags			{ IN_MAIN = 1, IN_SECN = 2 }; // Main, or Secondary, or both?
	void			put_in_MAIN()		{ x |= IN_MAIN; }
	void			put_in_SECN()		{ x |= IN_SECN; }
	bool			test_in_MAIN() const	{ return x & IN_MAIN; }
	bool			test_in_SECN() const	{ return x & IN_SECN; }
	typedef      bool (TriState :: *test_function_type) () const;
};

long double		empty_one_cluster(
					const int cluster_to_empty,
					Score & sc
				) {
	long double delta_score = 0.0L;
	const std :: tr1 :: unordered_set<int64_t> the_edges = sc.state.get_comms().at(cluster_to_empty).get_my_edges();
	For(edge, the_edges) {
		delta_score += sc.remove_edge(*edge, cluster_to_empty);
	}
	return delta_score;

}
long double		set_up_launch_state(
					const int main_cluster,
					const int secondary_cluster,
					const vector<int64_t>	& edges_in_a_random_order,
					Score & sc
				) {
	// Initially, ensure none of the edges are in either of the two communities
	For(edge, edges_in_a_random_order) {
		const int64_t edge_id = *edge;
		assert( sc.state.get_comms().at(main_cluster).get_my_edges().count(edge_id) == 0);
		assert( sc.state.get_comms().at(secondary_cluster).get_my_edges().count(edge_id) == 0);
	}


	const double alpha[3] = {1.0, 1.0, 1.0};
	double theta[3];
	gsl_ran_dirichlet(r, 3, alpha, theta);
	// for(int i=0;i<3;++i) { PP2(i, theta[i]); } cout << endl;

	long double delta_score = 0.0L;
	For(edge, edges_in_a_random_order) {
		unsigned int n[3];
		gsl_ran_multinomial(r, 3, 1, theta, n);
		// PP3(n[0],n[1],n[2]);
		if(n[0]) {
			delta_score += sc.add_edge(*edge, main_cluster);
		}
		if(n[1]) {
			delta_score += sc.add_edge(*edge, secondary_cluster);
		}
		if(n[2]) {
			delta_score += sc.add_edge(*edge, main_cluster);
			delta_score += sc.add_edge(*edge, secondary_cluster);
		}
	}
	// cout << endl;


	// Just one Gibbs sweep for now.  I think Jain&Neal showed that five sweeps was good for them though
	// Let's force the first two nodes into different clusters MUST BE TESTED THOROUGHLY
	for(int jain_neal_repeats=0; jain_neal_repeats<5; ++jain_neal_repeats) {
		For(edge, edges_in_a_random_order) {
			delta_score += gibbsUpdateJustTwoComms(*edge, sc, main_cluster, secondary_cluster).first;
		}
	}
	return delta_score;
}
long double		set_up_launch_state_POLICY_WITHIN_SPLIT_MERGE(
					const int main_cluster,
					const int secondary_cluster,
					const vector<int64_t>	& edges_in_a_random_order,
					Score & sc
				) {
	if(gsl_ran_bernoulli(r, 0.5)) {
		return set_up_launch_state(main_cluster, secondary_cluster, edges_in_a_random_order, sc);
	} else {
		return 0.0L;
	}
}

typedef std :: tr1 :: unordered_map< int64_t , TriState > Original_state_of_these_edges_T;
Original_state_of_these_edges_T remember_the_state_of_these_edges(const State &st, const int main_cluster, const int secondary_cluster) {
	Original_state_of_these_edges_T original_state_of_these_edges;
	For(main_edge, st.get_comms().at(main_cluster     ).get_my_edges()) { original_state_of_these_edges[ *main_edge ] . put_in_MAIN(); }
	For(secn_edge, st.get_comms().at(secondary_cluster).get_my_edges()) { original_state_of_these_edges[ *secn_edge ] . put_in_SECN(); }
	return original_state_of_these_edges;
}
vector<int64_t> those_edges_in_a_random_order(const Original_state_of_these_edges_T & original_state) {
	vector<int64_t> edges_in_a_random_order; // each edge to appear just once
	For(edge_with_state, original_state) {
		edges_in_a_random_order.push_back(edge_with_state->first);
	}
	assert(original_state.size() == edges_in_a_random_order.size());
	random_shuffle(edges_in_a_random_order.begin(), edges_in_a_random_order.end());
	return edges_in_a_random_order;
}

pair<long double, long double> from_launch_state_to_forced_proposal(
					const vector<int64_t>			& edges_in_a_random_order,
					int64_t					  main_cluster,
					int64_t					  secondary_cluster,
					      Original_state_of_these_edges_T	& original_state_of_these_edges,
					Score					& sc
		) {
				long double delta_score = 0.0L;
				long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES = 0.0L;
				For(edge, edges_in_a_random_order) {
					pair<int,int> possibly_force(0,0);
					if(original_state_of_these_edges [*edge].test_in_MAIN())
						possibly_force.first = 1;
					if(original_state_of_these_edges [*edge].test_in_SECN())
						possibly_force.second = 1;
					pair<long double, long double> ab = gibbsUpdateJustTwoComms(*edge, sc, main_cluster, secondary_cluster, possibly_force);
					log2_product_of_accepted_probabilities_FOR_ALL_EDGES += ab.second;
					delta_score += ab.first;
				}
				return make_pair(delta_score, log2_product_of_accepted_probabilities_FOR_ALL_EDGES);
}
pair<long double, long double> from_launch_state_to_UNforced_proposal(
					const vector<int64_t>			& edges_in_a_random_order,
					int64_t					  main_cluster,
					int64_t					  secondary_cluster,
					Score					& sc
		) {
				long double delta_score = 0.0L;
				long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES = 0.0L;
				For(edge, edges_in_a_random_order) {
					pair<long double, long double> ab = gibbsUpdateJustTwoComms(*edge, sc, main_cluster, secondary_cluster);
					log2_product_of_accepted_probabilities_FOR_ALL_EDGES += ab.second;
					delta_score += ab.first;
				}
				return make_pair(delta_score, log2_product_of_accepted_probabilities_FOR_ALL_EDGES);
}

long double		move_from_secondary_into_main(Score &sc, const int64_t main_cluster, const int64_t secondary_cluster, const vector<int64_t> & edges_in_a_random_order) {
		long double delta_score = 0.0L;
		// Remove every edge from the secondary_cluster, if it wasn't already
		delta_score += empty_one_cluster(secondary_cluster, sc);

		// Put every edge into the main_cluster, if it wasn't already
		For(edge, edges_in_a_random_order) {
			delta_score += sc.add_edge_if_not_already(*edge, main_cluster);
		}

		assert( sc.state.get_comms().at(secondary_cluster).get_my_edges().size() == 0);
		return delta_score;
}

pair<int, int> two_distinct_clusters(const int64_t K) {
	assert(K>=2);
	const int main_cluster = K * gsl_rng_uniform(r);
	int secondary_cluster;
	do { secondary_cluster = K * gsl_rng_uniform(r); }
	while (secondary_cluster == main_cluster);
	assert(secondary_cluster != main_cluster);
	return make_pair(main_cluster, secondary_cluster);
}

long double		M3(Score &sc) { // Does not change K
	if(sc.state.get_K() < 2)
		return 0.0L;
	// - Select two clusters at random
	// - Remember their current state
	// - Randomize the order of the edges

	// - Empty both of them, and Set up launch state
	// - Remember the launch state
	// - Do the "proposal", but with "forcing" of course

	// - Now, empty again and set up launch state again
	// - Reapply the launch state
	// - Do the  proposal , without forcing

	// - Calculate acceptance probability, and proceed as usual

	long double delta_score = 0.0L;
	const long double pre_score = delta_score;

	// - Select two clusters at random
	pair<int, int> two_clusters = two_distinct_clusters(sc.state.get_K());
	const int main_cluster = two_clusters.first;
	const int secondary_cluster = two_clusters.second;

	// - Remember their current state
	Original_state_of_these_edges_T original_state_of_these_edges = remember_the_state_of_these_edges(sc.state, main_cluster, secondary_cluster);

	// - Randomize the order of the edges
	const vector<int64_t> edges_in_a_random_order = those_edges_in_a_random_order(original_state_of_these_edges);

	// - Empty both of them, and Set up launch state
	delta_score += empty_one_cluster(main_cluster, sc);
	delta_score += empty_one_cluster(secondary_cluster, sc);

	/* Do not use a launch state here in M3. It actually makes it perform very badly.
	 * We need an easy-ish way to get to the bad (pre-) state.
	//delta_score += set_up_launch_state(main_cluster, secondary_cluster, edges_in_a_random_order      , sc);
	 */
	const long double score_at_launch_state = delta_score;
	Original_state_of_these_edges_T launch_state_of_these_edges = remember_the_state_of_these_edges(sc.state, main_cluster, secondary_cluster);

	// - Do the "proposal", but with "forcing" of course
	pair<long double, long double> result_forced = from_launch_state_to_forced_proposal(
			edges_in_a_random_order,
			main_cluster,
			secondary_cluster,
			original_state_of_these_edges,
			sc);
	const long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED = result_forced.second;
	delta_score += result_forced.first;

	// - Now, reapply the launch state
		delta_score += empty_one_cluster(main_cluster, sc);
		delta_score += empty_one_cluster(secondary_cluster, sc);
		For(edge_with_state, launch_state_of_these_edges) {
			const int64_t edge_id = edge_with_state->first;
			if(edge_with_state->second.test_in_MAIN())
				delta_score += sc.add_edge_if_not_already(edge_id, main_cluster);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, main_cluster);
			if(edge_with_state->second.test_in_SECN())
				delta_score += sc.add_edge_if_not_already(edge_id, secondary_cluster);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, secondary_cluster);

		}
		assertVERYCLOSE(delta_score, score_at_launch_state);

	// - Do the  proposal , without forcing
	pair<long double, long double> result_unforced = from_launch_state_to_UNforced_proposal(
			edges_in_a_random_order,
			main_cluster,
			secondary_cluster,
			sc);
	const long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED = result_unforced.second;
	delta_score += result_unforced.first;

	const long double post_score = delta_score;

	// - Calculate acceptance probability, and proceed as usual
	const long double acceptance_prob = post_score - pre_score
				- log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED
				+ log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED;
	//PP5(pre_score, post_score, log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED, log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED, acceptance_prob);

	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		// Accept this new split
		//PP(acceptance_prob);
		return delta_score;
	} else {
		// Rejected. Must return to the original state
		For(edge_with_state, original_state_of_these_edges) {
			const int64_t edge_id = edge_with_state->first;
			if(edge_with_state->second.test_in_MAIN())
				delta_score += sc.add_edge_if_not_already(edge_id, main_cluster);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, main_cluster);
			if(edge_with_state->second.test_in_SECN())
				delta_score += sc.add_edge_if_not_already(edge_id, secondary_cluster);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, secondary_cluster);

		}
		assertVERYCLOSE(delta_score, 0.0L);
		return delta_score;
	}
}

long double merge_these_two(Score &sc, const int main_cluster, const int secondary_cluster, const long double adjustment_to_acceptance = 0.0L) {
	assert(main_cluster >= 0);
	assert(secondary_cluster >= 0);
	assert(main_cluster != secondary_cluster);
	assert(secondary_cluster < sc.state.get_K());
	assert(main_cluster < sc.state.get_K());

	long double delta_score = 0.0L;

	// Remember their current state:
	Original_state_of_these_edges_T original_state_of_these_edges = remember_the_state_of_these_edges(sc.state, main_cluster, secondary_cluster);

	// Randomize the order of the edges:
	const vector<int64_t> edges_in_a_random_order = those_edges_in_a_random_order(original_state_of_these_edges);

	// Before emptying them, let's calculate the score of the merged version
	delta_score += move_from_secondary_into_main(sc, main_cluster, secondary_cluster, edges_in_a_random_order);
	assert( sc.state.get_comms().at(main_cluster).get_my_edges().size() == original_state_of_these_edges.size());
	const long double this_is_the_merged_score = delta_score + sc.what_would_change_if_I_deleted_an_empty_community();

	// Now, empty both of them:
	{
		delta_score += empty_one_cluster(main_cluster, sc);
		delta_score += empty_one_cluster(secondary_cluster, sc); // This is already empty, but it's no harm to be explicit

		assert(sc.state.get_comms().at(secondary_cluster).empty());
		assert(sc.state.get_comms().at(main_cluster     ).empty());
	}


	// Set up the launch state BUT ONLY 50% OF THE TIME
	delta_score += set_up_launch_state_POLICY_WITHIN_SPLIT_MERGE(main_cluster, secondary_cluster, edges_in_a_random_order      , sc);

	// Now finally ready for the "random" proposal
	pair<long double, long double> result = from_launch_state_to_forced_proposal(
			edges_in_a_random_order,
			main_cluster,
			secondary_cluster,
			original_state_of_these_edges,
			sc);
	const long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES = result.second;
	delta_score += result.first;

	// We're back at the start again. Unmerged!
	assertVERYCLOSE(delta_score, 0.0L);

	// Now, calculate the acceptance probability
	//   The cmf at the target (merging)        is this_is_the_merged_score {relatively speaking}
	//   The cmf at the source (split/original) is 0.0L                     {relatively speaking}
	//   The proposal probabilities are {up to proportionality}:
	//   		source -> target	0
	//   		target -> source	log2_product_of_accepted_probabilities_FOR_ALL_EDGES
	//   {The proposal probabilities associated with selecting a pair of clusters will cancel}
	// We accept the merge with acceptance probability exp2l{ this_is_the_merged_score + log2_product_of_accepted_probabilities_FOR_ALL_EDGES }
	const long double acceptance_prob = this_is_the_merged_score
		+ log2_product_of_accepted_probabilities_FOR_ALL_EDGES
		+ adjustment_to_acceptance;
	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		// Accept
		// This means I must merge them again
		{
			delta_score += move_from_secondary_into_main(sc, main_cluster, secondary_cluster, edges_in_a_random_order);

			assert( sc.state.get_comms().at(secondary_cluster).get_my_edges().size() == 0);
			assert( sc.state.get_comms().at(main_cluster).get_my_edges().size() == original_state_of_these_edges.size());

			// Now, to delete that empty cluster
			sc.state.swap_cluster_to_the_end(secondary_cluster);
			delta_score += sc.delete_empty_cluster_from_the_end();
			assertVERYCLOSE(delta_score, this_is_the_merged_score);
		}
		return delta_score;
	} else {
		// Leave them unmerged
		assertVERYCLOSE(delta_score, 0.0L);
		return delta_score;
	}
}
long double		split(Score &sc, const bool adjust_for_shared_edge_proposal = false) {
	// - Select one cluster at random, main_cluster
	// - Select another cluster_id at random, secondary_cluster, to be the cluster id of the new cluster
	// - (They do *not* need to different from each other)
	// - Randomize the order of the edges
	// - Empty it, create the new (empty) secondary_cluster, and Set up launch state
	// - Do the proposal, (no need for "forcing" here)
	// - Calculate acceptance probability, and proceed as usual

	const int main_cluster      =  sc.state.get_K()    * gsl_rng_uniform(r);
	const int secondary_cluster_NOT_TO_BE_USED_UNTIL_THE_SPLIT_IS_ACCEPTED = (sc.state.get_K()+1) * gsl_rng_uniform(r);

	long double delta_score = 0.0L;
	/// cout << endl << "About to attempt a split\t";
	/// dump_all(sc.state);
	const long double this_is_the_merged_score = delta_score; // == 0.0;

	// Gather the list of edges.
	// Randomize the order of the edges:
	vector<int64_t> edges_in_a_random_order; // each edge to appear just once
	For(main_edge, sc.state.get_comms().at(main_cluster     ).get_my_edges()) {
		edges_in_a_random_order.push_back(*main_edge);
	}
	assert(sc.state.get_comms().at(main_cluster).get_my_edges().size() == edges_in_a_random_order.size());
	random_shuffle(edges_in_a_random_order.begin(), edges_in_a_random_order.end());

	// Now, empty it:
	{
		For(edge, edges_in_a_random_order) {
			delta_score += sc.remove_edge(*edge, main_cluster);
		}
		assert(sc.state.get_comms().at(main_cluster     ).empty());
	}
	// Create the new, empty, secondary_cluster
	const int temporary_secondary_cluster_on_the_end = sc.state.get_K();
	delta_score += sc.append_empty_cluster(); // initially, this is on the end.  We'll only swap it in if the split it accepted.
	assert(temporary_secondary_cluster_on_the_end < sc.state.get_K());


	// Set up the launch state BUT ONLY 50% OF THE TIME
	delta_score += set_up_launch_state_POLICY_WITHIN_SPLIT_MERGE(main_cluster, temporary_secondary_cluster_on_the_end, edges_in_a_random_order      , sc);

	// Now finally ready for the "random" proposal
	long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES = 0.0L;
	For(edge, edges_in_a_random_order) {
		pair<long double, long double> ab = gibbsUpdateJustTwoComms(*edge, sc, main_cluster, temporary_secondary_cluster_on_the_end);
		log2_product_of_accepted_probabilities_FOR_ALL_EDGES += ab.second;
		delta_score += ab.first;
	}

	long double adjustment_to_acceptance = 0.0L;
	if(adjust_for_shared_edge_proposal) {
		// Now, we have our split.  What's the probability of proposing the reverse edge?
		// For example, their must be a shared edge.
		const long double p_shared_edge = probability_of_selecting_these_two_comms(main_cluster, temporary_secondary_cluster_on_the_end, sc.state);
		//  If there's no shared edge, then this will be minus-infinity and it will work out fine below.
		const int bigK = sc.state.get_K();
		const int smallK = bigK - 1;
		adjustment_to_acceptance += p_shared_edge + log2l(bigK) + log2l(smallK);
	}

	// Store the cmf {relatively speaking} at this split state
	const long double this_is_the_split_score = delta_score;

	// Now, calculate the acceptance probability
	//   The cmf at the target (split)           is this_is_the_split_score {relatively speaking}
	//   The cmf at the source (merged/original) is 0.0L                    {relatively speaking}
	//   The proposal probabilities are {up to proportionality}:
	//   		source -> target	log2_product_of_accepted_probabilities_FOR_ALL_EDGES
	//   		target -> source	0
	//   {The proposal probabilities associated with selecting a pair of clusters will cancel}
	// We accept the split with acceptance probability exp2l{ this_is_the_split_score - log2_product_of_accepted_probabilities_FOR_ALL_EDGES }
	const long double acceptance_prob = this_is_the_split_score - log2_product_of_accepted_probabilities_FOR_ALL_EDGES + adjustment_to_acceptance;
	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		// Accept
		// The split is accepted.  Not much more to do, except to swap
		sc.state.swap_cluster_to_the_end(secondary_cluster_NOT_TO_BE_USED_UNTIL_THE_SPLIT_IS_ACCEPTED);
		return delta_score;
	} else {
		// Rejected.  Must merge them again.
		{
			// delta_score += empty_one_cluster(temporary_secondary_cluster_on_the_end, & TriState :: test_in_SECN, original_state_of_these_edges, sc);
			// Put every edge into the main_cluster, and not in the secondary, if they weren't already
			// cout << "About to remerge\t";
			// dump_all(sc.state);
			For(edge, edges_in_a_random_order) {
				if( sc.state.get_edge_to_set_of_comms().at(*edge).count(main_cluster) == 0)
					delta_score += sc.add_edge(*edge, main_cluster);
				if( sc.state.get_edge_to_set_of_comms().at(*edge).count(temporary_secondary_cluster_on_the_end) == 1)
					delta_score += sc.remove_edge(*edge, temporary_secondary_cluster_on_the_end);
			}
			assert( sc.state.get_comms().at(temporary_secondary_cluster_on_the_end).get_my_edges().size() == 0);
			assert( sc.state.get_comms().at(main_cluster).get_my_edges().size() == edges_in_a_random_order.size());

			// Now, to delete that empty cluster
			// cout << "remerged\t";
			// dump_all(sc.state);
			delta_score += sc.delete_empty_cluster_from_the_end();
			// cout << "deleted the temp\t";
			// dump_all(sc.state);
			assertVERYCLOSE(delta_score, this_is_the_merged_score);
		}
		assertVERYCLOSE(delta_score, 0.0L);
		return delta_score;
	}
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

pair<int, int> find_two_comms_that_share_an_edge(const State &st) {
	const int64_t random_edge = st.E * gsl_rng_uniform(r);
	const vector<int64_t> comms_at_this_edge (st.get_edge_to_set_of_comms().at(random_edge).begin(), st.get_edge_to_set_of_comms().at(random_edge).end());
	assert( ! comms_at_this_edge.empty() );
	assert( comms_at_this_edge.size() == st.get_edge_to_set_of_comms().at(random_edge).size() );
	const int num_comms_here = comms_at_this_edge.size();
	if( num_comms_here == 1 ) {
		return make_pair(-1,-1);
	}
	assert(num_comms_here >= 2);
	const int main_cluster = num_comms_here * gsl_rng_uniform(r);
	int secondary_cluster;
	do {
		secondary_cluster = num_comms_here * gsl_rng_uniform(r);
	} while(secondary_cluster == main_cluster);
	return make_pair(main_cluster, secondary_cluster);
}

struct most_negative_ {
	operator long double() {
		return - numeric_limits<long double> :: max();
	}
};
most_negative_ most_negative() { return most_negative_(); }

static long double probability_of_selecting_these_two_comms(const int main_cluster, const int secondary_cluster, const State &st) {
	// First, find out which edges are shared between these two clusters.
	vector<int64_t> shared_edges;
	For(main_edge, st.get_comms().at(main_cluster).get_my_edges()) {
		if( st.get_comms().at(secondary_cluster).get_my_edges().count(*main_edge)) {
			shared_edges.push_back(*main_edge);
		}
	}
	if( shared_edges.empty() )
		return most_negative(); // log2 of zero is effectively minus-infinity
	long double total_probability = 0.0L;
	For(shared_edge, shared_edges) {
		const int64_t num_comms_at_this_shared_edge = st.get_edge_to_set_of_comms().at(*shared_edge).size();
		assert(num_comms_at_this_shared_edge >= 2);
		total_probability += 1.0L / num_comms_at_this_shared_edge / (num_comms_at_this_shared_edge-1);
	}
	assert(isfinite(total_probability));
	return log2l(total_probability);
}

static long double		merge_two_random_clusters(Score &sc) {
	if(sc.state.get_K() < 2)
		return 0.0L;
	// - Select two clusters at random
	// - Remember their current state
	// - Randomize the order of the edges
	// - Before emptying them, let's calculate the score of the merged version
	// - Empty both of them, and Set up launch state
	// - Do the "proposal", but with "forcing" of course
	// - Calculate acceptance probability, and proceed as usual

	pair<int, int> two_clusters = two_distinct_clusters(sc.state.get_K());
	const int main_cluster = two_clusters.first;
	const int secondary_cluster = two_clusters.second;
	return merge_these_two(sc, main_cluster, secondary_cluster);
}
long double		split_or_merge(Score & sc) {
	if(gsl_ran_bernoulli(r, 0.5))
		return merge_two_random_clusters(sc);
	else
		return split(sc);
}
long double		split_or_merge_on_a_shared_edge(Score & sc) {
	// The differences are:
	// - When calculating acceptance, we put (p_shared_edge * K * K') into the mix
	// - We decide which two to merge based on a random edge
	if(gsl_ran_bernoulli(r, 0.5)) {
		const pair<int,int> two_comms = find_two_comms_that_share_an_edge(sc.state);
		if(two_comms.first == -1)
			return 0.0L;
		else {
			const long double p_shared_edge = probability_of_selecting_these_two_comms(two_comms.first, two_comms.second, sc.state);
			const int K = sc.state.get_K();
			const long double adjustment_to_acceptance = - p_shared_edge - log2l(K) - log2l(K-1);
			return merge_these_two(sc, two_comms.first, two_comms.second, adjustment_to_acceptance);
		}
	} else {
		return split(sc, true);
	}
}
