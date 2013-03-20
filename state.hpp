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

#include "network.hpp"

struct Community {
private:
	std :: tr1 :: unordered_set<int64_t> my_edges;
	std :: tr1 :: unordered_multiset<int64_t> my_nodes; // each node can appear multiple times
	int64_t unique_nodes_in_this_community;
public:
	Community() : unique_nodes_in_this_community(0) {}
	void add_edge(int64_t e, Net net) {
		const bool inserted = this->my_edges.insert(e).second;
		assert(inserted);
		const network :: EdgeSet :: Edge & edge = net->edge_set->edges.at(e);
		this->add_node(edge.left);
		if(edge.left!=edge.right)
			this->add_node(edge.right);
	}
	void add_node(int64_t n) {
		this->my_nodes.insert(n);
		const size_t how_many = this->my_nodes.count(n);
		assert(how_many>0);
		if(how_many == 1) { // we've just added this new unique node
			++ this->unique_nodes_in_this_community;
		}
	}
};

struct Communities {
	// This is just a vector<Community> that will grow, and maybe shrink, every now and then.
private:
	std :: vector<Community> comms;
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
