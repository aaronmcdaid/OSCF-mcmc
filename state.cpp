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
