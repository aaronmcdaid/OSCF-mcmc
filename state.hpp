/*
 * The OSCF state is based on the edges, not the nodes.
 * I should be able to do directed and undirected, but I'll focus on undirected for now.
 *
 * But I'm not sure yet whether to have overlapping *edges* or not;
 * I'll try to handle both. Allowing an edge in more than one community
 * is more correct, but a bit more awkward.
 *
 * - There shall be K communities.
 * - Each community shall know which edges are in it.
 * - Each community shall know which nodes are in it, so that we know its order().
 * -   That requires a map : (community,node) -> int , which knows how 'often' a given node is in the community
 * - We need to be able to 'swap' two communities
 *
 * As edges are added/removed to communities, the above structures will be updated.
 * But we still need to know, for a given edge, which community(ies) it is in. Either:
 *
 * - a map: edge_id -> community_id
 * OR
 * - a map: edge_id -> set(community_id)  (A non-empty set)
 *
 *
 * Finally, I will (sort of) allow "missing data" - some edges might be totally unassigned at certain times.
 *
 */
#include <vector>
#include <cassert>
#include "tr1/unordered_set"
#include <cmath>
#include <gsl/gsl_sf.h>

#include "network.hpp"

struct Community {
private:
	std :: tr1 :: unordered_set<int64_t> my_edges;
	std :: tr1 :: unordered_multiset<int64_t> my_nodes; // each node can appear multiple times
	int64_t num_unique_nodes_in_this_community;
public:
	Community() : num_unique_nodes_in_this_community(0) {}
	void add_edge(int64_t e, Net net) {
		const bool inserted = this->my_edges.insert(e).second;
		assert(inserted);
		const network :: EdgeSet :: Edge & edge = net->edge_set->edges.at(e);
		this->add_node(edge.left);
		assert(edge.left!=edge.right);
		// if there are self loops, don't forget to consider the final *two* factors in the four factors of f(m_k, s'_k)
		if(edge.left!=edge.right)
			this->add_node(edge.right);
	}
	void add_node(int64_t n) {
		this->my_nodes.insert(n);
		const size_t how_many = this->my_nodes.count(n);
		assert(how_many>0);
		if(how_many == 1) { // we've just added this new unique node
			++ this->num_unique_nodes_in_this_community;
		}
	}
	int64_t		get_num_edges()					const { return this->my_edges.size(); }
	int64_t		get_num_unique_nodes_in_this_community()	const { return this->num_unique_nodes_in_this_community; }
};

struct Communities {
	// This is just a vector<Community> that will grow, and maybe shrink, every now and then.
private:
	std :: vector<Community> comms; // the length of this vector is *not* necessarily equal to K
public:
	Community & at(const int64_t k) {
		if(k >= this->comms.size())
			this->comms.resize(k+1);
		return this->comms.at(k);
	}
};

struct State {
	Net net;
	const int64_t N;
	const int64_t E;

	int64_t K; // number of communities
	Communities comms;


	explicit State(Net net_);
};



static inline double LOG2GAMMA(const double x) {
        assert(x>0);
        return M_LOG2E * gsl_sf_lngamma(x);
}
static inline double LOG2GAMMA(const long double x) {
        return LOG2GAMMA(static_cast<double>(x));
}
static inline double LOG2FACT(const int x) {
        assert(x>0);
        return M_LOG2E * gsl_sf_lnfact(x);
}
static inline double LOG2BINOM(const int n, const int m) {
	assert(m<=n);
	assert(m>=0);
	return LOG2FACT(n) - LOG2FACT(m) - LOG2FACT(n-m);
}

struct Score { // every modification *should* go through here eventually, so as to track the score.
	// But for now, this is just a passive object that recalculates all the scores from scratch each time
	State &		state;
	explicit	Score(State & state_)	: state(state_) {
	}
	long double	score()			const {
							return this->prior_on_K();
	}
	long double	prior_on_K()		const { return -LOG2FACT(state.K); }
	long double	product_on_fs()		const {
							long double s = 0.0L;
							for(int k=0; k<state.K; ++k) {
								const Community & comm = state.comms.at(k);
								s += f(comm.get_num_edges(), comm.get_num_unique_nodes_in_this_community());
							}
							return s;
	}
	long double	f(const int64_t num_edges, const int64_t num_unique_nodes_in_this_community)	const {
							long double s_one_comm = 0.0L;
							for(int64_t sz = num_unique_nodes_in_this_community; sz<num_unique_nodes_in_this_community + 5 && sz<state.N; ++sz) {
								// four factors
								s_one_comm += - log2l(state.N+1);

								const int64_t num_pairs = sz * (sz-1) / 2; // will change if self-loops allowed
								s_one_comm += - log2l(1+num_pairs);

								{
									const int64_t s_prime = num_unique_nodes_in_this_community;
									const int64_t s = sz;
									s_one_comm   += LOG2BINOM(state.N - s_prime, s - s_prime)
									              - LOG2BINOM(state.N          , s          );
								}

								assert(num_edges >  0); // MUSTDO >=
								assert(num_edges <  num_pairs); // MUSTDO <=



								s_one_comm += LOG2BINOM(num_pairs, num_edges);
							}
							return s_one_comm;
	}
};
