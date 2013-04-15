#include "state.hpp"

#include<iostream>
#include<iomanip>
using namespace std;

#include "macros.hpp"

#include "format_flag_stack/format_flag_stack.hpp"
format_flag_stack :: FormatFlagStack stack;


		State :: State(Net net_) : net(net_), N(net_->N()), E(net_->E()), edge_to_set_of_comms(net_->E())
					   , frequencies_of_edge_occupancy(10000, 0)
					   , total_count_of_edge_assignments(0)
{
	this->K = 0;
	this->frequencies_of_edge_occupancy.at(0) = this->E;
	// initially, there are no communities, and no edges are in any community.
	// But we'll initialize by placing each each in a community of its own.
}
void		Community :: dump_me()			const	{
						cout
							<< this->num_unique_nodes_in_this_community
							<< ' ' << this->my_edges.size();
						{
							const int64_t possible_pairs = this->num_unique_nodes_in_this_community * (this->num_unique_nodes_in_this_community-1) / 2;
							assert((int64_t)this->my_edges.size() <= possible_pairs);
							const long double density = double(this->my_edges.size()) / double(possible_pairs);
							if(possible_pairs>0) {
								cout << ' ';
								cout << stack.push << fixed <<  setprecision(1);
								cout << 100.0L * density << "% ";
								cout << stack.pop;
							}
							else
								cout << ' ';
						}
}
long double	OneCommunitySummary :: density()			const {
			const int64_t possible_pairs = this->num_unique_nodes_in_this_community * (this->num_unique_nodes_in_this_community-1) / 2;
			if(possible_pairs == 0)
				return -1;
			return double(this->num_edges) / double(possible_pairs);
}
std :: vector<int64_t>	Community :: get_my_nodes_NO_COUNT() const {
		vector<int64_t> the_nodes;
		For(node_and_count, this->my_nodes) {
			the_nodes.push_back(node_and_count->first);
		}
		return the_nodes;
}
void			State :: swap_cluster_to_the_end(const int64_t cluster_id)	{
					assert(cluster_id < this->K);
					if(cluster_id + 1 == this->K)
						return; // nothing to swap.  Just silently return.
					const int64_t last_cluster_id = this->K - 1;
					assert(cluster_id != last_cluster_id);
					std :: tr1 :: unordered_set<int64_t> left_edges = this->comms.at(cluster_id)      .my_edges;
					std :: tr1 :: unordered_set<int64_t> last_edges = this->comms.at(last_cluster_id) .my_edges;

					// Don't forget that some edges might already be in both
					// For efficiency, I should ignore the intersection of those two sets
					For(left_edge, left_edges) {
						this->remove_edge( *left_edge, cluster_id);
					}
					For(last_edge, last_edges) {
						this->remove_edge( *last_edge, last_cluster_id);
					}
					For(left_edge, left_edges) {
						this->add_edge    ( *left_edge, last_cluster_id);
					}
					For(last_edge, last_edges) {
						this->add_edge    ( *last_edge, cluster_id);
					}
}
