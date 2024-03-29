#include "moves.hpp"

#include <vector>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
using namespace std;
using namespace lvalue_input;

#include "macros.hpp"

gsl_rng * r = NULL;
gsl_rng * rng() { return r; }
#define log2_one_plus_l(x) (M_LOG2E * log1pl(x))
/*extern*/ int64_t GLOBAL_constraint_min_K;

struct most_negative_ {
	operator long double() {
		return - numeric_limits<long double> :: max();
	}
};
most_negative_ most_negative() { return most_negative_(); }

static long double probability_of_selecting_these_two_comms(const int main_cluster, const int secondary_cluster, const State &st);

void			seed_the_random_number_generator(int seed) {
					assert( r == NULL );
					r = gsl_rng_alloc (gsl_rng_taus);
					gsl_rng_set(r, seed);

					srand(seed); // I *think* this will seed std::random_shuffle
}

static pair< pair<bool,bool>, long double>	bernoullis_not_all_failed_2(
							const long double p_k_1
							, const long double p_k_2
							, std :: pair<int,int> possibly_force = std :: make_pair(-1,-1)
							) {
		if(possibly_force.first != -1) {
			assert(possibly_force.first  >= 0 && possibly_force.first  < 2); // It's just zero or one
			assert(possibly_force.second >= 0 && possibly_force.second < 2); // It's just zero or one
		}
		// Four possibilities:
		//   *0 0 (but this isn't allowed!)
		//    0 1
		//    1 0
		//    1 1

		assert(p_k_1 > 0);
		assert(p_k_2 > 0);
		const long double prob_both_failing = (1.0L-p_k_1) * (1.0L-p_k_2);
		long double prob_at_least_one_success = 1.0L - prob_both_failing;
		if(prob_at_least_one_success < p_k_1) {
			assertVERYCLOSE(prob_at_least_one_success, p_k_1);
			prob_at_least_one_success = p_k_1;
		}
		assert(prob_at_least_one_success >= p_k_1); // BROKEN?

		// First, calculate the conditional probability of cluster 1 being assigned
		const long double cond_prob_A = p_k_1 / prob_at_least_one_success;
		unless(isfinite(cond_prob_A)) {
			PP(__LINE__, cond_prob_A, p_k_1, prob_at_least_one_success);
		}
		assert(isfinite(cond_prob_A));

		long double log2_product_of_accepted_probabilities = 0.0L;

		bool bA;
		if(possibly_force.first == -1) {
			bA = gsl_rng_uniform(r) < cond_prob_A;
		} else {
			bA = possibly_force.first;
		}
		//PP(bA, cond_prob_A, log2l(cond_prob_A), log2_one_plus_l(-cond_prob_A));
		log2_product_of_accepted_probabilities += bA ? log2l(cond_prob_A) : log2_one_plus_l(-cond_prob_A);
		unless(isfinite(log2_product_of_accepted_probabilities)) {
			if(!isfinite(log2_product_of_accepted_probabilities)) {
				log2_product_of_accepted_probabilities = most_negative();
				assertVERYCLOSE(cond_prob_A, 1);
				assert(isfinite(log2_product_of_accepted_probabilities));
			}
			assert(isfinite(log2_product_of_accepted_probabilities));
		}
		assert(isfinite(log2_product_of_accepted_probabilities));

		const long double cond_prob_B = bA ? p_k_2 : 1.0L;
		bool bB;
		if(possibly_force.second == -1) {
			bB = gsl_rng_uniform(r) < cond_prob_B;
		} else {
			bB = possibly_force.second;
		}
		log2_product_of_accepted_probabilities += bB ? log2l(cond_prob_B) : log2_one_plus_l(-cond_prob_B);
		unless(isfinite(log2_product_of_accepted_probabilities)) {
			if(isinf(log2_product_of_accepted_probabilities) == 1) {
				log2_product_of_accepted_probabilities = most_negative();
				assert(cond_prob_B == 1);
			}
		}
		assert(isfinite(log2_product_of_accepted_probabilities));

		assert(bA || bB);
		return make_pair( make_pair(bA,bB), log2_product_of_accepted_probabilities);
}

static pair< vector<bool>, long double>	bernoullis_not_all_failed(
							const vector<long double> &p_k
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

static long double		remove_edge_from_all_its_communities(int64_t e, Score &sc) {
		long double delta_in_here = 0.0L;
		while( ! sc.state.get_edge_to_set_of_comms().at(e).empty() ) {
			int64_t comm_id_to_remove = * sc.state.get_edge_to_set_of_comms().at(e).begin();
			delta_in_here += sc.remove_edge(e, comm_id_to_remove);
		}
		return delta_in_here;
}
static long double		remove_edge_from_one_community_if_present(int64_t e, Score &sc, const int64_t comm_id_to_remove) {
		if( sc.state.get_edge_to_set_of_comms().at(e).count(comm_id_to_remove) == 1 )
			return sc.remove_edge(e, comm_id_to_remove);
		else
			return 0.0L;
}

static long double		calculate_p_based_on_the_log_ratio(const long double delta_score_one_edge) {
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
static vector<int64_t>		findNearbyCommunities(in<State> state, const int64_t e) {
	// For this edge, find all the neighbouring edges,
	// i.e. edges which share one endpoint.  Ignoring this edge, of course.
	// Then return all the communities that are on those neighbouring edges.
	vector<int64_t> nearby_communities;
	Net net = state->net;
	const network :: EdgeSet :: Edge edge_details = net->edge_set->edges.at(e);
	in< vector<int> > left_juncs = net->i.at(edge_details.left).my_junctions;
	For(junc_id, *left_juncs) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(*junc_id);
		if(junc.edge_id != e) {
			For(comm_nearby, state->get_edge_to_set_of_comms().at(junc.edge_id)) {
				nearby_communities.push_back(*comm_nearby);
			}
		}
	}
	in< vector<int> > right_juncs = net->i.at(edge_details.right).my_junctions;
	For(junc_id, *right_juncs) {
	//For(junc_id, net->i.at(edge_details.right).my_junctions) {
		const network :: Junction junc = net->junctions->all_junctions_sorted.at(*junc_id);
		if(junc.edge_id != e) {
			For(comm_nearby, state->get_edge_to_set_of_comms().at(junc.edge_id)) {
				nearby_communities.push_back(*comm_nearby);
			}
		}
	}
	sort(nearby_communities.begin(), nearby_communities.end());
	nearby_communities.erase(
		unique(nearby_communities.begin(), nearby_communities.end())
		, nearby_communities.end()
	);
	return nearby_communities;
}
long double 		gibbsUpdateNearby(Score& sc, int64_t e) {
	// Does NOT change K
	//
	// 1. remove the edge from its *nearby* communities
	// 2. reassign to each in turn, calculating the delta-score in each case, giving the probability for that Bernoulli
	// 3. draw from the Bernoullis, but conditioning that it must be assigned to at least one community.

	// First, find all this edges *nearby* communities
	const vector<int64_t> nearby_communities = findNearbyCommunities(sc.state, e);

	const size_t K_nearby = nearby_communities.size();
	if(K_nearby == 0)
		return 0.0L;

	long double delta_score = 0.0L;

	// Then, remove this edge from all those nearby_communities
	For(nearby_comm, nearby_communities) {
		delta_score += sc.set(e, *nearby_comm, false);
	}

	const network :: EdgeSet :: Edge edge_details = sc.state.net->edge_set->edges.at(e);

	// For each of the nearby commmunities, how does the corresponding f change with the (re)addition of this edge?
	vector<long double> p_k_nearby;
	{
		For(nearby_comm, nearby_communities) {
			// add_edge will do:
			// - increase num_edges by 1
			// - increase num_nodes by 0, 1 or 2
			const long double expected_delta_score = sc.if_this_edge_is_added(*nearby_comm, edge_details.left, edge_details.right);

			//const long double delta_score_one_edge = sc.add_edge(e, *nearby_comm);
			//assertEQ(num_nodes_post, that_comm->get_num_unique_nodes_in_this_community());
			//sc.state.remove_edge(e, *nearby_comm);
			//assertEQ(num_nodes_pre, that_comm->get_num_unique_nodes_in_this_community());

			//assertEQ(expected_delta_score, delta_score_one_edge);
			p_k_nearby.push_back(calculate_p_based_on_the_log_ratio(expected_delta_score));
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
			const long double delta_score_one_edge = sc.if_this_edge_is_added(k, e);
			p_k.at(justIterateOverTwo) = calculate_p_based_on_the_log_ratio(delta_score_one_edge);
		}

	}

	// Assign the new values

	const pair< pair<bool,bool>,long double > new_values_for_this_edge = bernoullis_not_all_failed_2(p_k.at(0),p_k.at(1), possibly_force);
	{
		if(new_values_for_this_edge.first.first)
				delta_in_gibbs += sc.add_edge(e, main_cluster);
		if(new_values_for_this_edge.first.second)
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
	while(!sc.state.get_comms().at(cluster_to_empty).get_my_edges().empty()) {
		delta_score += sc.remove_edge(*sc.state.get_comms().at(cluster_to_empty).get_my_edges().begin(), cluster_to_empty);
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

static
long double merge_these_two(Score &sc, const int main_cluster, const int secondary_cluster, const long double adjustment_to_acceptance) {
	assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
	if(sc.state.get_K() == GLOBAL_constraint_min_K) {
		return 0.0L;  // we can't merge them if they're already at the minimum allowed value of K
	}
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
		assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
		return delta_score;
	} else {
		// Leave them unmerged
		assertVERYCLOSE(delta_score, 0.0L);
		assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
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
					assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
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
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
							return delta_score;
						} else {
							// Reject
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
							return 0.0L;
						}
					} else {
						// cout << "Attempt removal" << endl;
						// Attempt to Remove an empty cluster

						// First, select a cluster at random to be our target.
						// If it's empty, just bail out immediately
						assert(sc.state.get_K() > 0);
						assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
						if(sc.state.get_K() == GLOBAL_constraint_min_K) {
							// Not allowed to remove, must Reject
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
							return 0.0L;
						}
						const int64_t target_cluster_id = sc.state.get_K() * gsl_rng_uniform(r);
						if(!sc.state.get_comms().at(target_cluster_id).empty()) {
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
							return 0.0L;
						}

						const long double hypothetical_delta_score = sc.what_would_change_if_I_deleted_an_empty_community();

						if( log2l(gsl_rng_uniform(r)) < hypothetical_delta_score ) {
							// Accept. Do the deletion
							sc.state.swap_cluster_to_the_end(target_cluster_id);
							const long double delta_score = sc.delete_empty_cluster_from_the_end();
							assert(delta_score > 0);
							assertEQ(delta_score, hypothetical_delta_score);
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
							return delta_score;
						} else {
							// Reject
							assert(1==2); // It should always Accept this proposal, due to the non-increasing prior on K
							assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
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
	assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
	if(sc.state.get_K() == GLOBAL_constraint_min_K) {
		return 0.0L;
	}
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
	return merge_these_two(sc, main_cluster, secondary_cluster, 0.0L);
}
long double		split_or_merge(Score & sc) {
	if(gsl_ran_bernoulli(r, 0.5))
		return merge_two_random_clusters(sc);
	else
		return split(sc);
}
vector<int> my_edges_in_a_random_order(State &st, const int k) {
	vector<int> edges_in_a_random_order; // each edge to appear just once
	For(main_edge, st.get_comms().at(k     ).get_my_edges()) {
		edges_in_a_random_order.push_back(*main_edge);
	}
	assert(st.get_comms().at(k).get_my_edges().size() == edges_in_a_random_order.size());
	random_shuffle(edges_in_a_random_order.begin(), edges_in_a_random_order.end());
	return edges_in_a_random_order;
}
void expand_seed(const int comm_source, const int comm_new, const int seed_edge, const vector<int> &E, Net net, long double &delta_score, Score &sc) {
	assert(sc.state.get_one_community_summary(comm_new).num_edges == 0);
	assert(sc.state.get_one_community_summary(comm_source).num_edges == (int64_t)E.size());
	// For the purposes of this, we *only* consider the edges in E. Only
	// those edges contribute towards the 'score' associated with a given frontier node
	//
	// Add nodes to the pair of nodes in {seed_edge}
	//
	// Maintain a frontier, a set of nodes that are connected to the growing community
	// This means a mapping from the frontier-node-id to the number of connections it has
	// .. and a priority_queue to find the most connected one each time
	//
	//
	//
	const set<int> E_set( E.begin(), E.end() );
	assert(E_set.size() == E.size());
	PP(E_set.size());

	map<int, size_t> how_many_frontier_edges_I_have;
	priority_queue< pair<long double, int> > noisy_queue; // with some randomness to tie break
	set<int> nodes_in_growing_comm;
	set<int> edges_in_growing_comm; // Not sure how I'll use this

	auto move_node_into_growing_comm = [&](int node) -> vector<int> {
		bool was_inserted = nodes_in_growing_comm.insert(node).second;
		assert(was_inserted);
		bool was_erased = how_many_frontier_edges_I_have.erase(node);
		assert(was_erased);
		// Iterate over the edges coming out of 'node'
		in< vector<int> > my_juncs = net->i.at(node).my_junctions;

		vector<int> edges_I_moved_into_the_comm;

		for(auto const junc_id : *my_juncs) {
			const network :: Junction junc = net->junctions->all_junctions_sorted.at(junc_id);
			assert(junc.this_node_id == node);
			auto const neigh_edge = junc.edge_id;
			if(E_set.count(neigh_edge)) {
				// OK, the far node's score can be increased ...
				// ... unless of course it is already fully in the expanding community
				const int far_node = junc.far_node_id;
				if(nodes_in_growing_comm.count(far_node) == 1) {
					const bool was_inserted = edges_in_growing_comm.insert(neigh_edge).second;
					assert(was_inserted);
					delta_score += sc.add_edge(neigh_edge, comm_new);
					delta_score += sc.remove_edge(neigh_edge, comm_source);
					edges_I_moved_into_the_comm.push_back(neigh_edge);
				} else {
					how_many_frontier_edges_I_have[far_node] ++;
					noisy_queue.push( { how_many_frontier_edges_I_have[far_node] + gsl_rng_uniform(r) * 0.001 , far_node });
				}
			}
		}
		return edges_I_moved_into_the_comm;
	};

	const network:: EdgeSet:: Edge & seed_Edge = net->edge_set->edges.at(seed_edge);
	//cout << "seed.left" << endl;
	how_many_frontier_edges_I_have[seed_Edge.left] = 0;
	move_node_into_growing_comm(seed_Edge.left);

	//cout << "seed.right" << endl;
	how_many_frontier_edges_I_have[seed_Edge.right] = 1;
	move_node_into_growing_comm(seed_Edge.right);

	//PP2(nodes_in_growing_comm.size(), edges_in_growing_comm.size());

	auto identify_next_node_to_add = [&]() -> int {
		//PP(noisy_queue.size());
		while(!noisy_queue.empty()) {
			auto const x = noisy_queue.top();
			noisy_queue.pop();
			if(how_many_frontier_edges_I_have.count(x.second)) {
				size_t const how_many = how_many_frontier_edges_I_have[x.second];
				//PP2(how_many, x.first);
				assert(how_many == floor(x.first));
				return x.second;
			}
		}
		return -1; // Nothing left to add
	};

	while(1) {
		const long double delta_score_before_this_node = delta_score;
		int const top_candidate = identify_next_node_to_add();
		if(top_candidate == -1)
			break;
		auto const edges_I_moved_into_the_comm = move_node_into_growing_comm(top_candidate);
		//const int current_size = nodes_in_growing_comm.size();
		//PP(current_size, edges_in_growing_comm.size(), current_size*(current_size-1)/2, delta_score);

		if(delta_score_before_this_node > 0 && delta_score < delta_score_before_this_node) {
			// should undo this addition and get out
			for(auto const edge_to_undo : edges_I_moved_into_the_comm) {
				delta_score += sc.add_edge(edge_to_undo, comm_source);
				delta_score += sc.remove_edge(edge_to_undo, comm_new);
			}
			PP(delta_score_before_this_node, delta_score);
			assertVERYCLOSE(delta_score_before_this_node, delta_score);
			break;
		}
		if(delta_score < 0 && delta_score < delta_score_before_this_node && nodes_in_growing_comm.size() > 5)
			break;

	}
	assert(sc.state.get_one_community_summary(comm_new).num_edges + sc.state.get_one_community_summary(comm_source).num_edges == (int64_t)E.size());
}
pair<int, long double> select_cluster_at_random_weighted_by_edge(State &st) {
	int total_y_kij = 0;
	for (auto const & cluster : st.get_comms()) {
		total_y_kij += cluster.get_num_edges();
	}
	assert( total_y_kij >= st.net->E() );
	int64_t random_latent_edge_offset = gsl_rng_uniform_int(r, total_y_kij);

	size_t k=0;
	while(true) {
		const int edges_in_k = st.get_one_community_summary(k).num_edges;
		if(edges_in_k > random_latent_edge_offset) {
			return { k, (edges_in_k+.0L) / total_y_kij };
		}
		++k;
		random_latent_edge_offset -= edges_in_k;
	}
}
long double             split_or_merge_by_seed_expansion(Score &sc) {
	// Attempt a merge or split
	// Record current state (very simple if we're just proposing to split)
	// Empty the community(ies)
	// Prepare launch state, via seed expansion
	// Either:
	//    Do proposal randomly (if splitting)
	//    Calculate reverse proposal probability (if merging)

	long double delta_score = 0.0L;
	const long double p_attempt_split = 1.0L / sc.state.get_K(); // if K=1, then it will not attempt a merge, of course
	cout << endl;
	PP(sc.state.get_K());
	auto calculate_prop_prob_from_seed = [&sc](
				const vector<int64_t> & launch_k
				, const int k
				, const vector<int64_t> & launch_l
				, const int l
				, const double (&alpha)[3]
				) -> long double {
			size_t count_asis = 0;
			size_t count_other= 0;
			size_t count_both = 0;
			auto do_one_edge_here = [&](int64_t const e, const int comm_asis, const int comm_other) {
				bool const still_in_asis = sc.state.test_edge(e, comm_asis);
				bool const also_in_other = sc.state.test_edge(e, comm_other);
				assert(still_in_asis || also_in_other); // must now be in one or the other, or both.
				if(still_in_asis && also_in_other)
					++ count_both;
				else
				if(still_in_asis)
					++ count_asis;
				else
				if(also_in_other)
					++ count_other;
				else
					assert(1==2); // shouldn't get here
			};
			assert(k < sc.state.get_K());
			assert(l < sc.state.get_K());
			for(auto const e : launch_k) {
				do_one_edge_here(e, k, l);
			}
			for(auto const e : launch_l) {
				do_one_edge_here(e, l, k);
			}
			PP(count_asis, count_other, count_both);
			//assert(verify_asis == count_asis);
			//assert(verify_other== count_other);
			//assert(verify_both == count_both);
			// Now everything is in one of three 'clusters', and we can use the collapsed formula
			//
			const long double l2_pp =
				LOG2GAMMA( alpha[0]+alpha[1]+alpha[2] )
				-LOG2GAMMA( alpha[0]+alpha[1]+alpha[2]
						+count_asis+count_other+count_both )
				+LOG2GAMMA( count_asis  + alpha[0] )
				+LOG2GAMMA( count_other + alpha[1] )
				+LOG2GAMMA( count_both  + alpha[2] )
				-LOG2GAMMA( alpha[0] )
				-LOG2GAMMA( alpha[1] )
				-LOG2GAMMA( alpha[2] )
				;
			PP(l2_pp);
			return l2_pp;
	}; // end lambda for the proposal prob, which will be useful in both split and merge
	if(gsl_ran_bernoulli(r, p_attempt_split)) {
		// Multiple parts to the proposal probability:
		// 0. Deciding to split, with probability 1/K
		// 1. Selecting a cluster at random, ACCORDING TO HOW MANY EDGES IT HAS
		// -. (append empty cluster, "q+1" )
		// -. (prepare launch state)
		// -. (make proposal)
		// 2. calculate proposal probability of second stage
		// 3. select a cluster to be swapped with "q+1"
		// that gives us three numbers for the proposal
		//
		// the reverse proposal probability is pretty simple, with *two* componnents:
		// 0. Deciding to merge, with probability (1- 1/K)
		// 1. just the probability of selecting two *distinct* clusters to merge
		auto const k_and_prob = select_cluster_at_random_weighted_by_edge(sc.state);
		PP(k_and_prob.first, k_and_prob.second);

		const int k = k_and_prob.first;
		const int qplus1 = sc.state.get_K();
		vector<int> E = my_edges_in_a_random_order(sc.state, k); // I don't really need them to be random, but what the hell
		if(E.size() < 2)
			return 0.0L;
		delta_score += sc.append_empty_cluster(); // initially, this is on the end.  We'll only swap it in if the split it accepted.
		assert(!E.empty());
		const size_t seed_edge = E.front();
		expand_seed(k, qplus1, seed_edge, E, sc.state.net, delta_score, sc);
		// Every edge is now in k or qplus1, but not both
		// This is the launch state, which we must now record
		const auto &    launch_state_qplus1_ = sc.state.get_comms().at(qplus1).get_my_edges();
		vector<int64_t> launch_state_qplus1  (launch_state_qplus1_.begin(), launch_state_qplus1_.end() );
		const auto &    launch_state_k_      = sc.state.get_comms().at(k).get_my_edges();
		vector<int64_t> launch_state_k       (launch_state_k_.begin(), launch_state_k_.end() );
		assert(launch_state_k.size() + launch_state_qplus1.size() == E.size());
		// Next, assign to one of the THREE possible states
		// first, remove them
		assert(int(E.size()) == sc.state.get_one_community_summary(k).num_edges + sc.state.get_one_community_summary(qplus1).num_edges);
		//auto empty_one_cluster = [&](const int one_cluster) -> void {
			//while(!sc.state.get_comms().at(qplus1).get_my_edges().empty()) { delta_score += sc.remove_edge( *sc.state.get_comms().at(qplus1).get_my_edges().begin(), qplus1); }
		//};
		delta_score += empty_one_cluster(k, sc);
		delta_score += empty_one_cluster(qplus1, sc);
		assert(0 == sc.state.get_one_community_summary(k).num_edges + sc.state.get_one_community_summary(qplus1).num_edges);

		const double alpha[3] = { staticcast(double,E.size()), 0.9, 0.1}; // as-is, other, both
		double theta[3];
		gsl_ran_dirichlet(r, 3, alpha, theta);
		PP(theta[0], theta[1], theta[2]);
		size_t verify_asis = 0;
		size_t verify_other= 0;
		size_t verify_both = 0;
		auto assign_one_edge_randomly = [&](const int edge_id, const int asIs, const int other, Score &sc) -> void {
			unsigned int n[3];
			gsl_ran_multinomial(r, 3, 1, theta, n);
			assert( n[0] + n[1] + n[2] == 1);
			if(n[0]) { delta_score += sc.add_edge( edge_id, asIs ); ++verify_asis; }  // as-is
			if(n[1]) { delta_score += sc.add_edge( edge_id, other ); ++verify_other; } // other
			if(n[2]) {                                  // both
				delta_score += sc.add_edge( edge_id, other );
				delta_score += sc.add_edge( edge_id, asIs );
				++verify_both;
			}
		};
		for(auto const e : launch_state_k      )  assign_one_edge_randomly(e, k, qplus1, sc);
		for(auto const e : launch_state_qplus1 )  assign_one_edge_randomly(e, qplus1, k, sc);
		//PP(verify_asis, verify_other, verify_both);
		assert(int(E.size()) <= sc.state.get_one_community_summary(k).num_edges + sc.state.get_one_community_summary(qplus1).num_edges);

		const long double l2_pp_post_launch = calculate_prop_prob_from_seed(launch_state_k, k, launch_state_qplus1, qplus1, alpha);
		PP(delta_score);
		const int new_size = sc.state.get_K();
		const long double l2_pp_all_forward = log2l(p_attempt_split)+ log2l(k_and_prob.second)+ l2_pp_post_launch+ log2l(1.0/new_size);
		PP(                                   log2l(p_attempt_split), log2l(k_and_prob.second), l2_pp_post_launch, log2l(1.0/new_size), l2_pp_all_forward);
		cout << "reverse pp:" << endl;
		const long double l2_pp_all_reverse = log2l( 1.0L - 1.0L/new_size )+ -log2l(new_size)-log2l(new_size-1);
		PP(                                   log2l( 1.0L - 1.0L/new_size ), -log2l(new_size)-log2l(new_size-1), l2_pp_all_reverse );

		const long double l2_acpt_prob = delta_score - l2_pp_all_forward + l2_pp_all_reverse;
		PP(l2_acpt_prob);
		if( log2l(gsl_rng_uniform(r)) < l2_acpt_prob) {
			cout << "  (seed-expand) SPLIT ACCEPT"
				<< endl;
			const int swap_with_me = gsl_rng_uniform_int(r, sc.state.get_K());
			sc.state.swap_cluster_to_the_end(swap_with_me);
			return delta_score;
		} else {
			// OK, but everything back the way it was
			PP(__LINE__, delta_score);
			for(const auto e : E) {
				if(!sc.state.test_edge(e, k))
					delta_score += sc.add_edge(e, k);
				if( sc.state.test_edge(e, qplus1))
					delta_score += sc.remove_edge(e, qplus1);
			}
			PP(__LINE__, delta_score);
			assert(int(E.size()) == sc.state.get_one_community_summary(k).num_edges);
			assert(  0           == sc.state.get_one_community_summary(qplus1).num_edges);
			PP(delta_score);
			delta_score += sc.delete_empty_cluster_from_the_end();
			PP(delta_score);
			assertVERYCLOSE(delta_score, 0);
			return delta_score;
		}
	}  else {
		assert(sc.state.get_K() > 1); // see 'p_attempt_split'

		// swap a cluster, randomly, to the end
		// DON'T FORGET TO UNDO THIS IF WE REJECT
		//
		// record the current state
		// empty qplus1, set k instead
		// set up the launch state
		// calculate the FORCED proposal probability
		// ...

		const int swap_with_me = gsl_rng_uniform_int(r, sc.state.get_K());
		sc.state.swap_cluster_to_the_end(swap_with_me);
		// how, we have our "expelled" cluster at the end
		const int k = gsl_rng_uniform_int(r, sc.state.get_K()-1);
		const int qplus1 = sc.state.get_K()-1;
		assert(k < qplus1);

		// record the current state
		const auto &    initial_state_qplus1_ = sc.state.get_comms().at(qplus1).get_my_edges();
		vector<int64_t> initial_state_qplus1  (initial_state_qplus1_.begin(), initial_state_qplus1_.end() );
		const auto &    initial_state_k_      = sc.state.get_comms().at(k).get_my_edges();
		vector<int64_t> initial_state_k       (initial_state_k_.begin(), initial_state_k_.end() );

		set<int> E_set( initial_state_k.begin(), initial_state_k.end() );
		assert(E_set.size() == initial_state_k.size());
		for(const auto e : initial_state_qplus1) {
			E_set.insert(e);
		}
		if(E_set.size() < 2) {
			sc.state.swap_cluster_to_the_end(swap_with_me);
			assertVERYCLOSE(delta_score, 0.0L);
			return delta_score; // == 0.0L
		}
		PP(E_set.size(), initial_state_qplus1.size(), initial_state_k.size());
		assert(E_set.size() <= initial_state_qplus1.size() + initial_state_k.size());
		vector<int> E (E_set.begin(), E_set.end());
		assert(E.size() == E_set.size());
		random_shuffle(E.begin(), E.end());

		// empty qplus1, put everything in k instead
		for(const auto edge_in_q_plus1 : initial_state_qplus1) {
			delta_score += sc.remove_edge(edge_in_q_plus1, qplus1);
			if(!sc.state.test_edge(edge_in_q_plus1, k))
				delta_score += sc.add_edge(edge_in_q_plus1, k);
		}
		const size_t seed_edge = E.front();
		expand_seed(k, qplus1, seed_edge, E, sc.state.net, delta_score, sc);

		// launch state now established
		// This is the launch state, which we must now record
		const auto &    launch_state_qplus1_ = sc.state.get_comms().at(qplus1).get_my_edges();
		vector<int64_t> launch_state_qplus1  (launch_state_qplus1_.begin(), launch_state_qplus1_.end() );
		const auto &    launch_state_k_      = sc.state.get_comms().at(k).get_my_edges();
		vector<int64_t> launch_state_k       (launch_state_k_.begin(), launch_state_k_.end() );
		assert(launch_state_k.size() + launch_state_qplus1.size() == E.size());
		assert(int(E.size()) == sc.state.get_one_community_summary(k).num_edges + sc.state.get_one_community_summary(qplus1).num_edges);

		// OK, go back to the 'initial' state!
		delta_score += empty_one_cluster(qplus1, sc);
		delta_score += empty_one_cluster(k, sc);
		for(const auto e : initial_state_k)      { delta_score += sc.add_edge(e, k); }
		for(const auto e : initial_state_qplus1) { delta_score += sc.add_edge(e, qplus1); }
		assertVERYCLOSE(delta_score, 0.0L);

		const double alpha[3] = { staticcast(double,E.size()), 0.9, 0.1}; // as-is, other, both
		const long double l2_pp_post_launch = calculate_prop_prob_from_seed(launch_state_k, k, launch_state_qplus1, qplus1, alpha);
		PP(l2_pp_post_launch);

		// Must force a merge now, to calculate the delta_score there, and the k_and_prob_second
		delta_score += empty_one_cluster(qplus1, sc);
		for(const auto e : E) {
			if(!sc.state.test_edge(e, k))
				delta_score += sc.add_edge(e, k);
		}
		assert((int64_t)E.size() == sc.state.get_one_community_summary(k).num_edges);

		delta_score += sc.delete_empty_cluster_from_the_end();
		long double k_and_prob_second;
		{
			size_t total_y_kij = 0;
			for(const auto & comm : sc.state.get_comms()) {
				total_y_kij += comm.get_num_edges();
			}
			k_and_prob_second = (sc.state.get_one_community_summary(k).num_edges+.0L) / total_y_kij;
		}

		// Now, accept or reject?
		// (currently, in the merged state)
		const int orig_larger_k = sc.state.get_K()+1;
		const int new_small_k = sc.state.get_K();
		cout << "reverse pp:" << endl;
		const long double l2_pp_all_reverse = log2l(1.0L/new_small_k)+ log2l(k_and_prob_second)+ l2_pp_post_launch+ log2l(1.0/orig_larger_k);
		PP(                                   log2l(1.0L/new_small_k), log2l(k_and_prob_second), l2_pp_post_launch, log2l(1.0/orig_larger_k), l2_pp_all_reverse);
		cout << "forward pp:" << endl;
		const long double l2_pp_all_forward = log2l( 1.0L - 1.0L/orig_larger_k )+ -log2l(orig_larger_k)-log2l(orig_larger_k-1);
		PP(                                   log2l( 1.0L - 1.0L/orig_larger_k ), -log2l(orig_larger_k)-log2l(orig_larger_k-1), l2_pp_all_forward );

		const long double l2_acpt_prob = delta_score - l2_pp_all_forward + l2_pp_all_reverse;
		PP(l2_acpt_prob);
		if( log2l(gsl_rng_uniform(r)) < l2_acpt_prob) {
			cout << "  (seed-expand) MERGE ACCEPT"
				<< endl;
			return delta_score;
		} else {
			// MERGE reject
			assert(qplus1 == sc.state.get_K());
			delta_score += sc.append_empty_cluster();
			// put the initial_state in play
			delta_score += empty_one_cluster(k, sc);
			for(const auto e : initial_state_k)
				delta_score += sc.add_edge(e, k);
			for(const auto e : initial_state_qplus1)
				delta_score += sc.add_edge(e, qplus1);
			// then swap again
			sc.state.swap_cluster_to_the_end(swap_with_me);
			assertVERYCLOSE(delta_score, 0.0L);
			return delta_score; // == 0.0L
		}
	}
}
long double		split_or_merge_on_a_shared_edge(Score & sc) {
	assert(sc.state.get_K() >= GLOBAL_constraint_min_K);
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
