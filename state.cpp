#include "state.hpp"

#include "macros.hpp"

		State :: State(Net net_) : net(net_), N(net_->N()), E(net_->E()), edge_to_set_of_comms(net_->E()) {
	this->K = 0;
	// initially, there are no communities, and no edges are in any community.
	// But we'll initialize by placing each each in a community of its own.
	for(int64_t e = 0; e < E; ++e) {
		const int64_t new_cluster_id = this->append_empty_cluster();
		assert(new_cluster_id == e);
		this->add_edge(e, e);
	}
}
void		Community :: dump_me()			const	{
						PP3(this->num_unique_nodes_in_this_community, this->my_edges.size(), this->my_nodes.size());
}
void			State :: swap_cluster_to_the_end(const int64_t cluster_id)	{
					assert(cluster_id < this->K);
					if(cluster_id + 1 == this->K)
						return; // nothing to swap.  Just silently return.
					const int64_t last_cluster_id = this->K - 1;
					assert(cluster_id != last_cluster_id);
					std :: tr1 :: unordered_set<int64_t> left_edges = this->comms.at(cluster_id)      .my_edges;
					std :: tr1 :: unordered_set<int64_t> last_edges = this->comms.at(last_cluster_id) .my_edges;

					For(left_edge, left_edges) {
						this->remove_edge( *left_edge, cluster_id);
						this->add_edge    ( *left_edge, last_cluster_id);
					}
					For(last_edge, last_edges) {
						this->remove_edge( *last_edge, last_cluster_id);
						this->add_edge    ( *last_edge, cluster_id);
					}
}
