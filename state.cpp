#include "state.hpp"
State :: State(Net net_) : net(net_), E(net_->E()), N(net_->N()) {
	this->K = 0;
	// initially, there are no communities, and no edges are in any community.
	// But we'll initialize by placing each each in a community of its own.
	for(int64_t e = 0; e < E; ++e) {
		this->comms.at(e).add_edge(e, net);
	}
}
