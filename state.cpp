#include "state.hpp"

State :: State(Net net_) : net(net_), E(net_->E()), N(net_->N()), edge_to_set_of_comms(net_->E()) {
	this->K = 0;
	// initially, there are no communities, and no edges are in any community.
	// But we'll initialize by placing each each in a community of its own.
	for(int64_t e = 0; e < E; ++e) {
		const int64_t new_cluster_id = this->append_empty_cluster();
		assert(new_cluster_id == e);
		this->add_edge(e, e);
	}
}
