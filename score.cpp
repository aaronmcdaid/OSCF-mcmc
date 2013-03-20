#include "score.hpp"

#include "macros.hpp"
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
							long double s_one_comm = 0.0L;
							for(int64_t sz = num_unique_nodes_in_this_community; sz<num_unique_nodes_in_this_community + 5 && sz<state.N; ++sz) {
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

								PP4(num_edges, num_unique_nodes_in_this_community, sz, score_one_comm_one_sz);

								s_one_comm += exp2l(score_one_comm_one_sz);
							}
							PP(log2l(s_one_comm));
							std :: cout << std :: endl;
							return log2l(s_one_comm);
}
