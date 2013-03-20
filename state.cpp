#include "state.hpp"

#include "macros.hpp"
State :: State(Net net_) : net(net_), E(net_->E()), N(net_->N()) {
	this->K = 0;
	// initially, there are no communities, and no edges are in any community.
	// But we'll initialize by placing each each in a community of its own.
	for(int64_t e = 0; e < E; ++e) {
		const int64_t new_cluster_id = this->append_empty_cluster();
		assert(new_cluster_id == e);
		this->add_edge(e, e);
	}
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
