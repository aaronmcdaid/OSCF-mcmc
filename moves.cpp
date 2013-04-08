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
using namespace lvalue_input;

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
static	std :: tr1 :: unordered_set<int64_t>		findNearbyCommunities(in<State> state, const int64_t e) {
	// For this edge, find all the neighbouring edges,
	// i.e. edges which share one endpoint.  Ignoring this edge, of course.
	// Then return all the communities that are on those neighbouring edges.
	std :: tr1 :: unordered_set<int64_t> nearby_communities; 
	Net net = state->net;
	const network :: EdgeSet :: Edge edge_details = net->edge_set->edges.at(e);
	For(junc_id, net->i.at(edge_details.left).my_junctions) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(*junc_id);
		if(junc.edge_id != e) {
			For(comm_nearby, state->get_edge_to_set_of_comms().at(junc.edge_id)) {
				nearby_communities.insert(*comm_nearby);
			}
		}
	}
	For(junc_id, net->i.at(edge_details.right).my_junctions) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(*junc_id);
		if(junc.edge_id != e) {
			For(comm_nearby, state->get_edge_to_set_of_comms().at(junc.edge_id)) {
				nearby_communities.insert(*comm_nearby);
			}
		}
	}
	return nearby_communities;
}
long double 		gibbsUpdateNearby(Score& sc, int64_t e) {
	// Does NOT change K
	//
	// 1. remove the edge from its *nearby* communities
	// 2. reassign to each in turn, calculating the delta-score in each case, giving the probability for that Bernoulli
	// 3. draw from the Bernoullis, but conditioning that it must be assigned to at least one community.

	// First, find all this edges *nearby* communities
	std :: tr1 :: unordered_set<int64_t> nearby_communities_ = findNearbyCommunities(sc.state, e);
	vector<int64_t> nearby_communities ( nearby_communities_.begin(), nearby_communities_.end());
	assert(nearby_communities.size() == nearby_communities_.size());

	const size_t K_nearby = nearby_communities.size();
	if(K_nearby == 0)
		return 0.0L;

	long double delta_score = 0.0L;

	// Then, remove this edge from all those nearby_communities
	For(nearby_comm, nearby_communities) {
		delta_score += sc.set(e, *nearby_comm, false);
	}


	// For each of the nearby commmunities, how does the corresponding f change with the (re)addition of this edge?
	vector<long double> p_k_nearby;
	{
		For(nearby_comm, nearby_communities) {
			const long double delta_score_one_edge = sc.add_edge(e, *nearby_comm);
			sc.state.remove_edge(e, *nearby_comm);
			p_k_nearby.push_back(calculate_p_based_on_the_log_ratio(delta_score_one_edge));
		}
		assert(p_k_nearby.size() == K_nearby);
	}

	if(!sc.state.get_edge_to_set_of_comms().at(e).empty()) {
		// This is easy, just do independent draws
		for(size_t k_nearby = 0; k_nearby < K_nearby; ++k_nearby) {
			const long double p_of_being_in_this_comm = p_k_nearby.at(k_nearby);
			assert(isfinite(p_of_being_in_this_comm));
			delta_score += sc.set(e, nearby_communities.at(k_nearby), gsl_rng_uniform(r) < p_of_being_in_this_comm);
		}
	} else {
		// Must use the conditional assignment
		const pair< vector<bool>,long double > new_values_for_this_edge = bernoullis_not_all_failed(p_k_nearby);
		for(size_t k_nearby = 0; k_nearby < K_nearby; ++k_nearby) {
			delta_score += sc.set(e, nearby_communities.at(k_nearby), new_values_for_this_edge.first.at(k_nearby));
		}
	}
	return delta_score;
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
	//PP3(theta[0], theta[1], theta[2]);

	long double delta_score = 0.0L;
	For(edge, edges_in_a_random_order) {
		unsigned int n[3];
		gsl_ran_multinomial(r, 3, 1, theta, n);
		assert( n[0] + n[1] + n[2] == 1);
		//PP3(n[0],n[1],n[2]);
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


	// I think Jain&Neal showed that five sweeps was good for them though
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
pair<int, int> two_distinct_clusters(const State &st) {
	return two_distinct_clusters(st.get_K());
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
	assert(isfinite(adjustment_to_acceptance));

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
						const long double hypothetical_delta_score = sc.what_would_change_if_I_added_an_empty_community();
						/// const long double post = sc.score();
						assert(hypothetical_delta_score < 0);
						/// assert(VERYCLOSE(post - pre, delta_score));
						if( log2l(gsl_rng_uniform(r)) < hypothetical_delta_score ) {
							// Accept
							const long double delta_score = sc.append_empty_cluster();
							assert(hypothetical_delta_score == delta_score);
							// .. but let's move it to a random location
							const int64_t target_cluster_id = sc.state.get_K() * gsl_rng_uniform(r);
							sc.state.swap_cluster_to_the_end(target_cluster_id);
							return delta_score;
						} else {
							// Reject
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

						const long double hypothetical_delta_score = sc.what_would_change_if_I_deleted_an_empty_community();

						if( log2l(gsl_rng_uniform(r)) < hypothetical_delta_score ) {
							// Accept. Do the deletion
							sc.state.swap_cluster_to_the_end(target_cluster_id);
							const long double delta_score = sc.delete_empty_cluster_from_the_end();
							assert(delta_score > 0);
							assertEQ(delta_score, hypothetical_delta_score);
							return delta_score;
						} else {
							// Reject
							assert(1==2); // It should always Accept this proposal, due to the non-increasing prior on K
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
	const int main_cluster_offset = num_comms_here * gsl_rng_uniform(r);
	int secondary_cluster_offset;
	do {
		secondary_cluster_offset = num_comms_here * gsl_rng_uniform(r);
	} while(secondary_cluster_offset == main_cluster_offset);
	return make_pair(comms_at_this_edge.at(main_cluster_offset), comms_at_this_edge.at(secondary_cluster_offset));
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
		total_probability += (1.0L / num_comms_at_this_shared_edge) / double(num_comms_at_this_shared_edge-1);
	}
	total_probability /= st.E;
	assert(isfinite(total_probability));
	assert(total_probability > 0.0L);
	assert(total_probability < 1.0L); // it should be less than (or close to) 0.5, I think
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
			assert(p_shared_edge != most_negative());
			const int K = sc.state.get_K();
			const long double adjustment_to_acceptance = - p_shared_edge - log2l(K) - log2l(K-1);
			return merge_these_two(sc, two_comms.first, two_comms.second, adjustment_to_acceptance);
		}
	} else {
		return split(sc, true);
	}
}

long double swap_this_node_wrt_two_clusters(Score &sc, const int64_t n, const int clusterA, const int clusterB, const vector<bool> &swap_or_not) {
	Net net = sc.state.net;
	long double delta_score = 0.0L;
	assert( swap_or_not.size() == net->i.at(n).my_junctions.size() );
	const int degree = swap_or_not.size();
	for(int d= 0; d < degree; ++d ) {
		const bool should_I_swap_this = swap_or_not.at(d);
		if(!should_I_swap_this)
			continue;
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(net->i.at(n).my_junctions.at(d));
		assert(n == junc.this_node_id);
		const int edge_id = junc.edge_id;

		const bool is_in_clusterA = sc.state.get_edge_to_set_of_comms().at(edge_id).count(clusterA);
		const bool is_in_clusterB = sc.state.get_edge_to_set_of_comms().at(edge_id).count(clusterB);
		if(is_in_clusterA && !is_in_clusterB) {
			delta_score += sc.remove_edge(edge_id, clusterA);
			delta_score += sc.add_edge(edge_id, clusterB);
		}
		if(is_in_clusterB && !is_in_clusterA) {
			delta_score += sc.remove_edge(edge_id, clusterB);
			delta_score += sc.add_edge(edge_id, clusterA);
		}
	}
	return delta_score;
}

long double one_node_CHEATING(Score &sc) {
	const int64_t N = sc.state.N;
	const int64_t random_node = gsl_rng_uniform(r) * N;

	pair<int, int> two_clusters = two_distinct_clusters(sc.state);
	const int clusterA = two_clusters.first;
	const int clusterB = two_clusters.second;

	const int degree = sc.state.net->i.at(random_node).total_degree();

	long double delta_score = 0.0L;
	const long double pre_score = delta_score;

	// randomize that one node's edge wrt both clusters
	// .. BUT ignoring edges that are BOTH outside the cluster
	// .. then perform a load of Gibbs updates
	// (I really should allow rejection afterwards, as per Jain+Neal

	const double alpha[3] = {1.0, 1.0, 1.0};
	double theta[3];
	gsl_ran_dirichlet(r, 3, alpha, theta);
	//PP3(theta[0], theta[1], theta[2]);

	Net net = sc.state.net;
	vector<int64_t> edges; // this is to be the set of all edges associated with this node, which are in at least one of the comms
	Original_state_of_these_edges_T original_state_of_these_edges;
	for(int d= 0; d < degree; ++d ) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(net->i.at(random_node).my_junctions.at(d));
		assert(random_node == junc.this_node_id);
		const int e = junc.edge_id;
		const bool is_in_clusterA = sc.state.get_edge_to_set_of_comms().at(e).count(clusterA);
		const bool is_in_clusterB = sc.state.get_edge_to_set_of_comms().at(e).count(clusterB);
		if(is_in_clusterA || is_in_clusterB)
			edges.push_back(e);
		if(is_in_clusterA) { original_state_of_these_edges[ e ] . put_in_MAIN(); }
		if(is_in_clusterB) { original_state_of_these_edges[ e ] . put_in_SECN(); }
	}
	assert(original_state_of_these_edges.size() == edges.size());
	if(edges.empty())
		return 0.0L;
	random_shuffle(edges.begin(), edges.end());
	assert((int)edges.size() <= degree);

	For(edge, edges) {
		const int e = *edge;

		assert(!sc.state.get_edge_to_set_of_comms().at(e).empty());
		delta_score += remove_edge_from_one_community_if_present(e, sc, clusterA);
		delta_score += remove_edge_from_one_community_if_present(e, sc, clusterB);

		unsigned int n[3];
		gsl_ran_multinomial(r, 3, 1, theta, n);
		assert( n[0] + n[1] + n[2] == 1);

		if(n[0]==1) { // both
			delta_score += sc.add_edge(e, clusterA);
			delta_score += sc.add_edge(e, clusterB);
		}
		if(n[1]==1) { /* just A */ delta_score += sc.add_edge(e, clusterA); }
		if(n[2]==1) { /* just B */ delta_score += sc.add_edge(e, clusterB); }
	}

	for(int rep=0; rep<20; ++rep) {
		For(edge, edges) {
			pair<long double, long double> ab = gibbsUpdateJustTwoComms(*edge, sc, clusterA, clusterB);
			delta_score += ab.first;
		}
	}

	const long double delta_score_at_launch_state = delta_score;

	// Store the launch state
	Original_state_of_these_edges_T launch_state_of_these_edges;
	For(edge, edges) {
		const int e = *edge;
		if(sc.state.get_edge_to_set_of_comms().at(e).count(clusterA)) { launch_state_of_these_edges[ e ] . put_in_MAIN(); }
		if(sc.state.get_edge_to_set_of_comms().at(e).count(clusterB)) { launch_state_of_these_edges[ e ] . put_in_SECN(); }
	}
	assert(launch_state_of_these_edges.size() == edges.size());

	// Force back to the original state, remembering the proposal probability
	pair<long double, long double> result_forced = from_launch_state_to_forced_proposal(
			edges,
			clusterA,
			clusterB,
			original_state_of_these_edges,
			sc);
	const long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED = result_forced.second;
	delta_score += result_forced.first;

	assertVERYCLOSE(delta_score, 0.0L);

	// Restore the launch state
	For(edge_with_state, launch_state_of_these_edges) {
		const int64_t edge_id = edge_with_state->first;
		if(edge_with_state->second.test_in_MAIN())
			delta_score += sc.add_edge_if_not_already(edge_id, clusterA);
		else
			delta_score += sc.remove_edge_if_not_already(edge_id, clusterA);
		if(edge_with_state->second.test_in_SECN())
			delta_score += sc.add_edge_if_not_already(edge_id, clusterB);
		else
			delta_score += sc.remove_edge_if_not_already(edge_id, clusterB);
	}
	assertVERYCLOSE(delta_score, delta_score_at_launch_state);

	pair<long double, long double> result_unforced = from_launch_state_to_UNforced_proposal(
			edges,
			clusterA,
			clusterB,
			sc);
	const long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED = result_unforced.second;
	delta_score += result_unforced.first;

	const long double post_score = delta_score;

	const long double acceptance_prob = post_score - pre_score
				- log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED
				+ log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED;

	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		// Accept. Do nothing
		return delta_score;
	} else {
		// Restore the original state
		For(edge_with_state, original_state_of_these_edges) {
			const int64_t edge_id = edge_with_state->first;
			if(edge_with_state->second.test_in_MAIN())
				delta_score += sc.add_edge_if_not_already(edge_id, clusterA);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, clusterA);
			if(edge_with_state->second.test_in_SECN())
				delta_score += sc.add_edge_if_not_already(edge_id, clusterB);
			else
				delta_score += sc.remove_edge_if_not_already(edge_id, clusterB);
		}
		assertVERYCLOSE(delta_score, pre_score);
		return delta_score;
	}

}

long double one_node_simple_update(Score &sc) {
	if(sc.state.get_K() == 1)
		return 0.0L;
	// Select a node at random,
	// and two distinct clusters at random
	// Swap the neighbouring edges from one cluster to another.
	const int64_t N = sc.state.N;
	const int64_t random_node = gsl_rng_uniform(r) * N;

	pair<int, int> two_clusters = two_distinct_clusters(sc.state);

	long double delta_score = 0.0L;

	// We won't swap all of them, just a randomly-selected subset of the edges
	const int degree = sc.state.net->i.at(random_node).total_degree();
	vector<bool> swap_or_not;
	const long double policy = gsl_rng_uniform(r);
	// PP(policy);
	for(int d = 0; d< degree; ++d) {
		swap_or_not.push_back( gsl_ran_bernoulli(r, policy) );
		// PP(swap_or_not.back());
	}
	assert(degree == (int)swap_or_not.size());

	delta_score += swap_this_node_wrt_two_clusters(sc, random_node, two_clusters.first, two_clusters.second, swap_or_not);

	const long double acceptance_prob = delta_score;
	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		// Accept, do nothing and return
		// PP2(random_node, sc.state.net->node_set->as_string(random_node));
	} else {
		// Reject, Undo
		delta_score += swap_this_node_wrt_two_clusters(sc, random_node, two_clusters.first, two_clusters.second, swap_or_not);
		assertVERYCLOSE(delta_score, 0.0L);
	}
	return delta_score;
}

int64_t get_edge( const pair<int64_t, int64_t> p ) { return p.first; }
int64_t get_opp ( const pair<int64_t, int64_t> p ) { return p.second; }

vector<bool> store_these_edges_wrt_this_comm(const vector< pair<int64_t,int64_t> > &edges, int comm, in<State> st) {
	vector<bool> is_in_here;
	For(edge, edges) {
		if(st->get_edge_to_set_of_comms().at( get_edge(*edge) ).count(comm))
			is_in_here.push_back(true);
		else
			is_in_here.push_back(false);
	}
	assert(is_in_here.size() == edges.size());
	return is_in_here;
}
vector<bool>	store_this_edge_wrt_these_comms(const int64_t edge, in< vector<int> > communities, in<State> st) {
	vector<bool> this_edge;
	For(comm, *communities) {
		const int comm_id = *comm;
		const bool b = st->get_edge_to_set_of_comms().at(edge).count(comm_id);
		this_edge.push_back(b);
	}
	assert(this_edge.size() == communities->size());
	return this_edge;
}

struct p_B_fourCases {
	const double p_B_opp_neither;
	const double p_B_opp_A;
	const double p_B_opp_B;
	const double p_B_opp_both;
	p_B_fourCases() : 
			p_B_opp_neither(gsl_rng_uniform(r))
			,p_B_opp_A(gsl_rng_uniform(r))
			,p_B_opp_B(gsl_rng_uniform(r))
			,p_B_opp_both(gsl_rng_uniform(r))
	{}
};
const pair<long double, long double> make_proposal_for_SIMPLEST(Score &sc
		, const int64_t edge, const int64_t opp
		, const int clusterA, const int clusterB
		, in<p_B_fourCases> p_B
		, const int possibly_force_into_B = 300 // 0 or 1 will signify to force into A or B
		) {
	assert(possibly_force_into_B == 300 || possibly_force_into_B == 0 || possibly_force_into_B == 1);
	const bool opp_is_in_A = sc.state.get_comms().at(clusterA).test_node(opp);
	const bool opp_is_in_B = sc.state.get_comms().at(clusterB).test_node(opp);
	double p_b = 2.0;
	if(opp_is_in_A && opp_is_in_B)   p_b = p_B->p_B_opp_both;
	if(opp_is_in_A && !opp_is_in_B)  p_b = p_B->p_B_opp_A;
	if(!opp_is_in_A && opp_is_in_B)  p_b = p_B->p_B_opp_B;
	if(!opp_is_in_A && !opp_is_in_B) p_b = p_B->p_B_opp_neither;
	assert(p_b < 1.1);

	long double delta_score = 0.0L;
	const bool b = possibly_force_into_B == 300 ? gsl_ran_bernoulli(r, p_b) : (possibly_force_into_B==1);
	delta_score += sc.set(edge, clusterA, !b);
	delta_score += sc.set(edge, clusterB, b);
	return make_pair(delta_score, b ? log2l(p_b) : log2l(1.0-p_b) );
}
long double one_node_SIMPLEST_update(Score &sc, const int64_t random_node) {
	// Select one node at random
	// Select two distinct communities at random
	// Find the edges of that node which are in EXACTLY one of those two communities.
	// There will be four types of opposite-endpoints:
	// 	Neither, A, B, Both
	// 	We need four probabilities, P(current -> A), for each of those possibilities.

	if(sc.state.get_K()==1)
		return 0.0L;

	Net net = sc.state.net;
	const int degree = sc.state.net->i.at(random_node).total_degree();

	pair<int,int> distinct = two_distinct_clusters(sc.state);
	const int clusterA = distinct.first;
	const int clusterB = distinct.second;

	vector< pair<int64_t,int64_t> > edges;
	vector<bool> original_state; // Is in B, not A
	for(int d= 0; d < degree; ++d ) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(net->i.at(random_node).my_junctions.at(d));
		assert(random_node == junc.this_node_id);
		const int e = junc.edge_id;
		const bool is_in_clusterA = sc.state.get_edge_to_set_of_comms().at(e).count(clusterA);
		const bool is_in_clusterB = sc.state.get_edge_to_set_of_comms().at(e).count(clusterB);
		if( is_in_clusterA != is_in_clusterB) { // in one, but not both
			const int opposite_end_node = junc.far_node_id;
			assert(random_node != opposite_end_node); // not sure how to handle self loops in this method.
			edges.push_back( make_pair(e, opposite_end_node) );
			original_state.push_back(is_in_clusterB);
		}
	}
	assert(original_state.size() == edges.size());
	if(edges.empty()) {
		return 0.0L;
	}

	const p_B_fourCases p_B;


	long double delta_score = 0.0L;

	long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED = 0.0L;
	for(size_t e=0; e<edges.size(); ++e) {
		const int64_t edge = get_edge( edges.at(e) );
		const int64_t opp = get_opp ( edges.at(e) );
		assert(opp != random_node);
		const pair<long double, long double> delta_and_propprob = make_proposal_for_SIMPLEST(sc, edge, opp, clusterA, clusterB, p_B, original_state.at(e));
		delta_score += delta_and_propprob.first;
		log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED += delta_and_propprob.second;
	}

	const long double score_forced = delta_score;
	assertVERYCLOSE(score_forced, 0.0L); // to verify we're back at the start

	long double log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED = 0.0L;
	for(size_t e=0; e<edges.size(); ++e) {
		const int64_t edge = get_edge( edges.at(e) );
		const int64_t opp = get_opp ( edges.at(e) );
		assert(opp != random_node);
		const pair<long double, long double> delta_and_propprob = make_proposal_for_SIMPLEST(sc, edge, opp, clusterA, clusterB, p_B);
		delta_score += delta_and_propprob.first;
		log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED += delta_and_propprob.second;
	}

	const long double score_unforced = delta_score;

	const long double acceptance_prob = score_unforced - score_forced - log2_product_of_accepted_probabilities_FOR_ALL_EDGES_UNFORCED
		+ log2_product_of_accepted_probabilities_FOR_ALL_EDGES_FORCED;

	if(log2l(gsl_rng_uniform(r)) < acceptance_prob) {
		//PP(acceptance_prob);
		return delta_score;
	} else {
		// Reject. Undo
		for(size_t e=0; e<edges.size(); ++e) {
			const int64_t edge = get_edge( edges.at(e) );
			bool orig_in_B = original_state.at(e);
			delta_score += sc.set( edge, clusterA, !orig_in_B);
			delta_score += sc.set( edge, clusterB,  orig_in_B);
		}
		assertVERYCLOSE(delta_score, 0.0L);
		return delta_score;
	}
}
long double		gibbs_one_comm_one_edge(Score & sc, const int64_t e) {
	const int K = sc.state.get_K();
	const int k = gsl_rng_uniform(r) * K;

	assert(!sc.state.get_edge_to_set_of_comms().at(e).empty());
	long double delta_score = 0.0L;
	delta_score += sc.set( e, k, false);

	// If that was the only edge, we just have to include it again :-(
	if(sc.state.get_edge_to_set_of_comms().at(e).empty()) {
		delta_score += sc.set( e, k, true);
		assertVERYCLOSE(0.0L, delta_score);
		assert(!sc.state.get_edge_to_set_of_comms().at(e).empty());
		return delta_score;
	}
	const long double extra_if_in = sc.set( e, k, true);
	delta_score += extra_if_in;

	const long double exponented = exp2l(extra_if_in);
	const long double prob_connnected = exponented / (1.0L + exponented);
	assert(isfinite(prob_connnected));

	if(gsl_rng_uniform(r) < prob_connnected) {
		// Connect
		// It's already connected, just return
		return delta_score;
	} else {
		// Disconnect
		delta_score += sc.set( e, k, false);
		return delta_score;
	}
}
