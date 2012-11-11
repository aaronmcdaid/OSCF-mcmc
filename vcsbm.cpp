const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"
gengetopt_args_info args_info; // a global variable! Sorry.
#include "macros.hpp"

#include "network.hpp"

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

	const network :: NodeSet * node_set = build_node_set_from_edge_list(edgeListFileName,
			args_info.stringIDs_flag
			?  network :: NODE_NAME_STRING
			:  network :: NODE_NAME_INT64
			);
	const network :: EdgeSet * edge_set = build_edge_set_from_edge_list(edgeListFileName,
			args_info.weighted_flag
			? network :: EdgeSet :: WEIGHT_INT
			: network :: EdgeSet :: WEIGHT_NONE
			, node_set
			);
	const int N = node_set -> N();
	const int E = edge_set->edges.size();
	PP2(N,E);

	// Finally, we create the list of Junctions - by doubling the list of Edges
	const network :: Junctions * junctions = build_junctions_set_from_edges(edge_set, args_info.directed_flag);
	PP(junctions->all_junctions_sorted.size());

#if 0 // This was the validation code, print the network out again to check it
	for(int i=0; i<10 && i<N; ++i) { // print the first ten nodes, to make sure the order is as expected.
		PP(node_set -> as_string(i));
	}
	for(size_t e = 0; e < edge_set->edges.size(); ++e) {
		string left_str = node_set->as_string(edge_set->edges.at(e).left);
		string right_str = node_set->as_string(edge_set->edges.at(e).right);
		cout << left_str << ' ' << right_str << endl;
	}
#endif
}
