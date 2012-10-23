#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <set>
#include "network.hpp"
#include "macros.hpp"

using namespace network;
using namespace std;

struct NodeSet_Impl : public NodeSet_I {
	virtual int get_NodeNameType() const {
		return this->node_name_type;
	}
	virtual int N() const {
		assert (this->node_name_type == NODE_NAME_INT64);
		if (this->node_name_type == NODE_NAME_INT64)
			return this->node_names_int64.size();
		return -1; // shouldn't get here
	}
	std :: string as_string(int) const {
		return "";
	}

	enum NodeNameType node_name_type;
	NodeSet_Impl(enum NodeNameType nnt) : node_name_type(nnt) {
	}

	// Only one of these two sets will be used.
	set<int64_t> node_names_int64;
	set<int64_t> node_names_string;

	void notify_of_node(int64_t node_name) {
		assert(this-> node_name_type == NODE_NAME_INT64);
		this->node_names_int64 . insert(node_name);
	}
};

NodeSet_I * network :: build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeSet_I :: NodeNameType node_name_type) {
	ifstream edgelist(edgeListFileName);
	if(!edgelist) {
		cerr << "Error: edge list file (" << edgeListFileName << ") not found. Exiting" << endl;
		exit(1);
	}
	NodeSet_Impl * nodes = new NodeSet_Impl(node_name_type);

	string line;
	int64_t line_num = 0;
	while (getline(edgelist, line)) {
		++ line_num;
		PP(line);
		int64_t left_node, right_node;
		istringstream fields(line);
		fields >> left_node;
		if(!fields) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
		fields >> right_node;
		if(!fields) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
		PP2(left_node, right_node);
		nodes -> notify_of_node(left_node);
		nodes -> notify_of_node(right_node);
	}
	return nodes;
}
