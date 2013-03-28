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
#ifndef STATE_HPP__
#define STATE_HPP__

#include <vector>
#include <cassert>
#include "tr1/unordered_set"
#include "tr1/unordered_map"
#include <gsl/gsl_sf.h>

#include "network.hpp"

struct State; // pre-declare in order that we can make it a friend of Community

struct OneCommunitySummary {
	int64_t num_edges;
	int64_t num_unique_nodes_in_this_community;
};

struct Community {
	friend struct State;
private:
	std :: tr1 :: unordered_set<int64_t> my_edges;
	std :: tr1 :: unordered_map<int64_t,int64_t> my_nodes; // each node can appear multiple times
	int64_t num_unique_nodes_in_this_community;
public:
	Community() : num_unique_nodes_in_this_community(0) {}
	void	dump_me()			const;
	const struct OneCommunitySummary get_one_community_summary() const {
		OneCommunitySummary ocs;
		ocs.num_edges = this->get_num_edges();
		ocs.num_unique_nodes_in_this_community = this->get_num_unique_nodes_in_this_community();
		return ocs;
	}
	std :: vector<int64_t> get_my_nodes_NO_COUNT() const;
private:
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
	void remove_edge(int64_t e, Net net) {
		const size_t wasErased = this->my_edges.erase(e);
		assert(wasErased == 1);
		const network :: EdgeSet :: Edge & edge = net->edge_set->edges.at(e);
		this->remove_node(edge.left);
		assert(edge.left!=edge.right);
		// if there are self loops, don't forget to consider the final *two* factors in the four factors of f(m_k, s'_k)
		if(edge.left!=edge.right) {
			this->remove_node(edge.right);
		}
	}
	void add_node(int64_t n) {
		this->my_nodes[n]++;
		const size_t how_many = this->my_nodes[n];
		assert(how_many>0);
		if(how_many == 1) { // we've just added this new unique node
			++ this->num_unique_nodes_in_this_community;
		}
	}
	void remove_node(int64_t n) {
		std :: tr1 :: unordered_map <int64_t,int64_t> :: iterator it = this->my_nodes.find(n);
		assert(it != this->my_nodes.end());
		assert(it->second>0);
		assert(it->first == n);
		it->second --;
		const size_t how_many = it->second;
		if(how_many == 0) { // we've just removed the last copy of this node
			-- this->num_unique_nodes_in_this_community;
			this->my_nodes.erase(it);
		}
	}
public:
	int64_t		get_num_edges()					const { return this->my_edges.size(); }
	int64_t		get_num_unique_nodes_in_this_community()	const { return this->num_unique_nodes_in_this_community; }
	bool		empty()						const {
										if(this->num_unique_nodes_in_this_community == 0) {
											assert(this->my_edges.empty());
											assert(this->my_nodes.empty());
											return true;
										} else {
											assert(this->num_unique_nodes_in_this_community > 0);
											assert(!this->my_edges.empty());
											assert(!this->my_nodes.empty());
											return false;
										}
	}
	const std :: tr1 :: unordered_set<int64_t>		& get_my_edges()	const { return this->my_edges; }
};

struct State {
	friend class Score; // Score is allowed to edit this
	friend void dump_all(const State &);
	Net net;
	const int64_t N;
	const int64_t E;

private:
	int64_t K; // number of communities
	std :: vector<Community> comms; // the length of this vector is *not* necessarily equal to K

	std :: vector< std :: tr1 :: unordered_set<int64_t> > edge_to_set_of_comms;

	std :: vector<size_t> frequencies_of_edge_occupancy;
	std :: int64_t total_count_of_edge_assignments;


public:
	explicit State(Net net_);
	int64_t								  get_K()			const { return this->K; }
	const std :: vector<Community>					& get_comms()			const { return this->comms; }
	const std :: vector< std :: tr1 :: unordered_set<int64_t> >	& get_edge_to_set_of_comms()	const { return this->edge_to_set_of_comms; }
	const OneCommunitySummary					  get_one_community_summary(const int k)	const {
		return this->comms.at(k).get_one_community_summary();
	}
	// I'll allow {add,remove}_edge to be public, as they return void
	// and don't imply that they calculate the delta-score.
	// If there are Score::{add,remove}_edge, they should calculate the delta-score
public:
	void		add_edge(int64_t e, int64_t comm_id)		{
										assert(comm_id < this->K);
										this->comms.at(comm_id).add_edge(e, this->net);
										-- this->frequencies_of_edge_occupancy.at( this->edge_to_set_of_comms.at(e).size() );
										bool wasInserted = this->edge_to_set_of_comms.at(e).insert(comm_id).second;
										assert(wasInserted);
										++ this->frequencies_of_edge_occupancy.at( this->edge_to_set_of_comms.at(e).size() );
										++ this->total_count_of_edge_assignments;
	}
	void		remove_edge(int64_t e, int64_t comm_id)		{
										assert(comm_id < this->K);
										this->comms.at(comm_id).remove_edge(e, this->net);
										-- this->frequencies_of_edge_occupancy.at( this->edge_to_set_of_comms.at(e).size() );
										size_t wasErased = this->edge_to_set_of_comms.at(e).erase(comm_id);
										assert(wasErased == 1);
										++ this->frequencies_of_edge_occupancy.at( this->edge_to_set_of_comms.at(e).size() );
										-- this->total_count_of_edge_assignments;
	}
private:
	int64_t		append_empty_cluster()				{
										assert(this->K == (int64_t) this->comms.size());
										const int64_t new_cluster_id = this->K;
										this->comms.push_back( Community() );
										assert(this->comms.at(new_cluster_id).get_num_edges() == 0);
										++ this->K;
										assert(this->K == (int64_t)this->comms.size());
										return new_cluster_id;
	}
	void		delete_empty_cluster_from_the_end()		{
										assert(this->K == (int64_t)this->comms.size());
										const int64_t cluster_id_to_delete = int64_t(this->K)-1;
										assert(cluster_id_to_delete >= 0);  // can't delete when there are no clusters left!
										assert(this->comms.at(cluster_id_to_delete).empty());
										-- this->K;
										this->comms.pop_back();
										assert(this->K == (int64_t)this->comms.size());
	}
public: // swap is public because it doesn't affect the score
	void		swap_cluster_to_the_end(const int64_t cluster_id);
	bool		every_edge_non_empty()		const {
				for(int64_t e = 0; e<net->E(); ++e) {
					if(this->edge_to_set_of_comms.at(e).empty())
						return false;
				}
				return true;
	}
};




#endif
