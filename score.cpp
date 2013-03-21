#include "score.hpp"

#include <cmath>

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
							return this->prior_on_K() + this->product_on_fs();
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
							long double s_one_comm_sans_baseline = 0.0L;
							long double baseline = NAN; // this 'baseline' technique should improve accuracy
							for(int64_t sz = num_unique_nodes_in_this_community; sz<num_unique_nodes_in_this_community + 10 && sz<state.N; ++sz) {
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
							}
							// PP2(baseline , log2l(s_one_comm_sans_baseline));
							// PP(baseline + log2l(s_one_comm_sans_baseline));
							// std :: cout << std :: endl;
							return baseline + log2l(s_one_comm_sans_baseline);
}
