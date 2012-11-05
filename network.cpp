#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <set>
#include <vector>
#include <unordered_map>
#include "network.hpp"
#include "macros.hpp"

using namespace network;
using namespace std;

struct NodeSet_Impl : public NodeSet_I {
	virtual int get_NodeNameType() const {
		return this->node_name_type;
	}
	virtual int N() const {
		assert(locked());
		return this -> N_;
	}
	std :: string as_string(int node_id) const {
		assert(node_id >= 0 && node_id < this->N());
		if (this->node_name_type == NODE_NAME_INT64) {
			ostringstream oss;
			oss << this->node_names_int64_sorted.at(node_id);
			return oss.str();
		} else
			return this->node_names_string_sorted.at(node_id);
	}

	enum NodeNameType node_name_type;
	NodeSet_Impl(enum NodeNameType nnt) : node_name_type(nnt) {
	}

	// Only one of these two sets will be used.
	set<int64_t> node_names_int64;
	set<string> node_names_string;

	void notify_of_node(int64_t node_name) {
		assert(!locked());
		assert(this-> node_name_type == NODE_NAME_INT64);
		this->node_names_int64 . insert(node_name);
	}
	void notify_of_node(string node_name) {
		assert(!locked());
		assert(this-> node_name_type == NODE_NAME_STRING);
		this->node_names_string . insert(node_name);
	}

	int N_;
	vector<int64_t> node_names_int64_sorted;
	vector<string> node_names_string_sorted;
	unordered_map<int64_t, int> map_int64_to_id;
	unordered_map<string, int> map_string_to_id;
	bool locked() const {
		return !(node_names_int64_sorted.empty() && node_names_string_sorted.empty());
	}
	void lock() {
		assert(!locked());
		if(this->node_name_type == NODE_NAME_INT64) {
			copy(node_names_int64.begin(), node_names_int64.end(), back_inserter(node_names_int64_sorted));
			this -> N_ = node_names_int64_sorted.size();
			for(int i=0; i<this->N(); ++i) {
				map_int64_to_id[ node_names_int64_sorted.at(i) ] = i;
			}
		} else {
			copy(node_names_string.begin(), node_names_string.end(), back_inserter(node_names_string_sorted));
			this -> N_ = node_names_string_sorted.size();
			for(int i=0; i<this->N(); ++i) {
				map_string_to_id[ node_names_string_sorted.at(i) ] = i;
			}
		}
		assert(locked());
		assert(node_names_int64_sorted.size() == map_int64_to_id.size());
		assert(node_names_string_sorted.size() == map_string_to_id.size());
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
		// PP(line);
		int64_t left_node_int64, right_node_int64;
		string left_node_string, right_node_string;
		istringstream fields(line);
		if(node_name_type == network :: NodeSet_I :: NODE_NAME_INT64)
			fields >> left_node_int64;
		else
			fields >> left_node_string;
		if(!fields) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
		if(node_name_type == network :: NodeSet_I :: NODE_NAME_INT64)
			fields >> right_node_int64;
		else
			fields >> right_node_string;
		if(!fields) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
		if(node_name_type == network :: NodeSet_I :: NODE_NAME_INT64) {
			// PP2(left_node_int64, right_node_int64);
			nodes -> notify_of_node(left_node_int64);
			nodes -> notify_of_node(right_node_int64);
		} else {
			// PP2(left_node_string, right_node_string);
			nodes -> notify_of_node(left_node_string);
			nodes -> notify_of_node(right_node_string);
		}
	}
	nodes->lock();
	return nodes;
}
