const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"
gengetopt_args_info args_info; // a global variable! Sorry.
#include "macros.hpp"

#include "network.hpp"

#include<algorithm>
#include<cassert>

using namespace std;

int main(int argc, char **argv) {
	// Parse the args - there should be exactly one arg, the edge list
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(1) ;
	if(args_info.git_version_flag) {
		PP(gitstatus);
		for (int i=0; i<argc; i++) {
			PP(argv[i]);
		}
	}
	if(args_info.inputs_num != 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	const char * edgeListFileName   = args_info.inputs[0];

	network :: NodeSet_I * node_set = build_node_set_from_edge_list(edgeListFileName,
			args_info.stringIDs_flag
			?  network :: NodeSet_I :: NODE_NAME_STRING
			:  network :: NodeSet_I :: NODE_NAME_INT64
			);
	network :: EdgeSet * edge_set = build_edge_set_from_edge_list(edgeListFileName,
			args_info.weighted_flag
			? network :: EdgeSet :: WEIGHT_INT
			: network :: EdgeSet :: WEIGHT_NONE
			, node_set
			);
	const int N = node_set -> N();
	const int E = edge_set->edges.size();
	PP2(N,E);
#if 0 // This was the validation code, print the network out again to check it
	for(size_t e = 0; e < edge_set->edges.size(); ++e) {
		string left_str = node_set->as_string(edge_set->edges.at(e).left);
		string right_str = node_set->as_string(edge_set->edges.at(e).right);
		cout << left_str << ' ' << right_str << endl;
	}
#endif
	// Finally, we create the list of Endpoints - by doubling the list of Edges
	struct EndPoint {
		int edge_id; // edge_id: a number between 0 and E
		int endpoint_type; // 1: I'm the source, and the other end of the sink
		                   // -1: Other way around
				   // 0: Undirected, so it doesn't matter.
		int this_node_id; // id of this node at this end of the edge.
		int far_node_id; // id of the node at the other end of this  edge.
		EndPoint(int edge_id_, int endpoint_type_, int this_node_id_, int far_node_id_)
			: edge_id(edge_id_), endpoint_type(endpoint_type_), this_node_id(this_node_id_), far_node_id(far_node_id_) {}
		bool operator <(const EndPoint &other) const {
			if(this->this_node_id < other.this_node_id) return true;
			if(this->this_node_id > other.this_node_id) return false;
			assert(this->this_node_id == other.this_node_id);
			if(this->far_node_id < other.far_node_id) return true;
			if(this->far_node_id > other.far_node_id) return false;
			assert(this->far_node_id == other.far_node_id);
			return false;
		}
	};
	vector<EndPoint> all_endpoints_sorted; // we'll sort this later
	for(int e=0; e<E; ++e) {
		const network :: EdgeSet :: Edge & edge = edge_set->edges.at(e);
		EndPoint l(e, args_info.directed_flag ? 1 : 0, edge.left, edge.right);
		EndPoint r(e, args_info.directed_flag ?-1 : 0, edge.right, edge.left);
		all_endpoints_sorted.push_back(l);
		all_endpoints_sorted.push_back(r);
	}
	sort(all_endpoints_sorted.begin(), all_endpoints_sorted.end());
}
