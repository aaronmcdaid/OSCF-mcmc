#include "score.hpp"

#include <cmath>
#include <limits>
using namespace std;

#include "macros.hpp"

static inline double LOG2GAMMA(const double x) {
        assert(x>0);
        return M_LOG2E * gsl_sf_lngamma(x);
}
static inline double LOG2GAMMA(const long double x) {
        return LOG2GAMMA(static_cast<double>(x));
}
static inline double LOG2FACT(const int x) {
        assert(x>=0);
        return M_LOG2E * gsl_sf_lnfact(x);
}
static inline double LOG2BINOM(const int n, const int m) {
	assert(m<=n);
	assert(m>=0);
	return LOG2FACT(n) - LOG2FACT(m) - LOG2FACT(n-m);
}

		Score :: Score(State & state_)	: state(state_) {
}
long double	Score :: score()			const {
							const long double total_score = this->prior_on_K() + this->product_on_fs();
							assert(isfinite(total_score));
							return total_score;
}
long double	Score :: prior_on_K()		const { return -LOG2FACT(state.K); }
long double	Score :: product_on_fs()		const {
							long double s = 0.0L;
							for(int k=0; k<state.K; ++k) {
								const Community & comm = state.comms.at(k);
								s += f(comm.get_num_edges(), comm.get_num_unique_nodes_in_this_community());
							}
							return s;
}
long double	Score :: f	(const int64_t num_edges, const int64_t num_unique_nodes_in_this_community)	const {
					// I should check for overflows here, and for conversion down to int32_t
							const Cache_T :: key_type key = make_pair(num_edges, num_unique_nodes_in_this_community);
							{
								// First, check the cache
								Cache_T :: iterator it = this->cache.find(key);
								if(it != this->cache.end()) {
									if(it->second.second < std :: numeric_limits< int64_t > :: max() )
										++ it->second.second; // increment the hit count
									assert(it->second.second > 0);
									return it->second.first;
								}
							}

							long double s_one_comm_sans_baseline = 0.0L;
							long double baseline = NAN; // this 'baseline' technique should improve accuracy
							assert(!isfinite(baseline));
							assert(num_unique_nodes_in_this_community <= state.N);
							for(int64_t sz = num_unique_nodes_in_this_community; sz<num_unique_nodes_in_this_community + 100 && sz<=state.N; ++sz) {
								long double score_one_comm_one_sz = 0.0L;
								const int64_t num_pairs = sz * (sz-1) / 2; // will change if self-loops allowed

								// four factors
								score_one_comm_one_sz += - log2l(state.N+1);

								{
									const int64_t s_prime = num_unique_nodes_in_this_community;
									const int64_t s = sz;
									assert(s_prime <= s);
									score_one_comm_one_sz   += LOG2BINOM(state.N - s_prime, s - s_prime)
									                          - LOG2BINOM(state.N          , s          );
								}

								score_one_comm_one_sz += - log2l(1+num_pairs);


								score_one_comm_one_sz -= LOG2BINOM(num_pairs, num_edges);

								assert(num_edges >= 0);
								assert(num_edges <= num_pairs);

								// PP4(num_edges, num_unique_nodes_in_this_community, sz, score_one_comm_one_sz);

								if(sz == num_unique_nodes_in_this_community) {
									// this is the first one, and score_one_comm_one_sz should therefore
									// be smaller than the upcoming ones
									baseline = score_one_comm_one_sz;
								}
								score_one_comm_one_sz -= baseline;

								s_one_comm_sans_baseline += exp2l(score_one_comm_one_sz);
								assert(isfinite(s_one_comm_sans_baseline));
							}
							assert(isfinite(baseline)); // That loop should succeed at least once
							// PP2(baseline , log2l(s_one_comm_sans_baseline));
							// PP(baseline + log2l(s_one_comm_sans_baseline));
							// std :: cout << std :: endl;
							const long double total = baseline + log2l(s_one_comm_sans_baseline);
							assert(isfinite(total));

							this->cache.insert( make_pair(key, make_pair(total,0) ) );
							return total;
}
long double	Score :: f	(const OneCommunitySummary ocs)	const {
	return this->f(ocs.num_edges, ocs.num_unique_nodes_in_this_community);
}

long double	Score :: add_edge(int64_t e, int64_t comm_id_to_add) {
			const OneCommunitySummary old_one_comm = this->state.get_one_community_summary(comm_id_to_add);
			this->state.add_edge(e, comm_id_to_add);
			const OneCommunitySummary new_one_comm = this->state.get_one_community_summary(comm_id_to_add);
			return this->f(new_one_comm) - this->f(old_one_comm);
}
long double	Score :: remove_edge(int64_t e, int64_t comm_id_to_remove) {
			const OneCommunitySummary old_one_comm = this->state.get_one_community_summary(comm_id_to_remove);
			this->state.remove_edge(e, comm_id_to_remove);
			const OneCommunitySummary new_one_comm = this->state.get_one_community_summary(comm_id_to_remove);
			return this->f(new_one_comm) - this->f(old_one_comm);
}

long double		Score :: append_empty_cluster()	 {
	const long double f_0_0 = this->f(0,0);
	const long double pre_prior = this->prior_on_K();
	this->state.append_empty_cluster();
	const long double post_prior = this->prior_on_K();
	assert(post_prior < pre_prior);
	return post_prior - pre_prior + f_0_0;
}
long double		Score :: delete_empty_cluster_from_the_end()	 {
	const long double verify = this->what_would_change_if_I_deleted_an_empty_community();
	const long double f_0_0 = this->f(0,0);
	const long double pre_prior = this->prior_on_K();
	this->state.delete_empty_cluster_from_the_end();
	const long double post_prior = this->prior_on_K();
	assert(post_prior > pre_prior);
	assert(verify == post_prior - pre_prior - f_0_0);
	return post_prior - pre_prior - f_0_0;
}
long double 		Score :: what_would_change_if_I_deleted_an_empty_community()				const {
	assert(this->state.K > 0);
	const long double f_0_0 = this->f(0,0);
	const long double pre_prior = this->prior_on_K();
	-- this->state.K;
	const long double post_prior = this->prior_on_K();
	assert(post_prior > pre_prior);
	++ this->state.K;
	return post_prior - pre_prior - f_0_0;
}
