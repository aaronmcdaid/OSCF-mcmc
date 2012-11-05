const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
#include "cmdline.h"
gengetopt_args_info args_info; // a global variable! Sorry.
#include "macros.hpp"

#include "network.hpp"

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
	PP(node_set -> N());
	for(int i = 0; i<node_set ->N(); ++i) {
		PP(node_set -> as_string(i));
	}
}
