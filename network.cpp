#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>
#include "network.hpp"
#include "macros.hpp"

using namespace network;
using namespace std;

struct EdgeListFileFormatException : public std :: exception {
};
struct EdgeListFileFormatException_ExpectedAnInt : public EdgeListFileFormatException {};

template<typename T>
struct NodeSet_ : public NodeSet {
	set<T> set_of_names;
	vector<T> vector_of_names;
	virtual void insert_string_version_of_name(string);
	virtual string as_string(int);
	virtual void finish_me();
	virtual int N() const;
	virtual NodeNameType get_NodeNameType() const;
};
	template<>
	void NodeSet_<string> :: insert_string_version_of_name(string s) { set_of_names.insert(s); }
	template<>
	void NodeSet_<int64_t> :: insert_string_version_of_name(string s) {
		int64_t name_as_int;
		istringstream iss(s);
		iss >> name_as_int;
		if(! ( iss && iss.eof() ) ) {
			throw EdgeListFileFormatException_ExpectedAnInt();
		}
		set_of_names.insert(name_as_int);
	}
	template<typename T>
	void NodeSet_<T> :: finish_me() {
		assert(!this->set_of_names.empty());
		copy(this->set_of_names.begin(), this->set_of_names.end(), back_inserter(vector_of_names));
		this->set_of_names.clear();
		assert(this->set_of_names.empty());
	}
	template<>
	string NodeSet_<string> :: as_string(int node_id) {
		assert(node_id >= 0);
		assert(node_id < (int)this->vector_of_names.size());
		return this->vector_of_names.at(node_id);
	}
	template<>
	string NodeSet_<int64_t> :: as_string(int node_id) {
		assert(node_id >= 0);
		assert(node_id < (int)this->vector_of_names.size());
		const int64_t node_name = this->vector_of_names.at(node_id);
		ostringstream oss;
		oss << node_name;
		return oss.str();
	}
	template<typename T>
	int NodeSet_<T> :: N() const {
		assert(this->set_of_names.empty());
		assert(!this->vector_of_names.empty());
		return this->vector_of_names.size();
	}
	template<> NodeNameType NodeSet_<int64_t> :: get_NodeNameType() const { return NODE_NAME_INT64; }
	template<> NodeNameType NodeSet_<string>  :: get_NodeNameType() const { return NODE_NAME_STRING; }
template struct NodeSet_<string>;
template struct NodeSet_<int64_t>;






struct NodeSet_Impl : public NodeSet_I {
	virtual NodeNameType get_NodeNameType() const {
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
	int64_t as_int64(int node_id) const {
		assert (this->node_name_type == NODE_NAME_INT64);
		return this->node_names_int64_sorted.at(node_id);
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
	bool locked() const {
		return !(node_names_int64_sorted.empty() && node_names_string_sorted.empty());
	}
	void lock() {
		assert(!locked());
		if(this->node_name_type == NODE_NAME_INT64) {
			copy(node_names_int64.begin(), node_names_int64.end(), back_inserter(node_names_int64_sorted));
			this -> N_ = node_names_int64_sorted.size();
		} else {
			copy(node_names_string.begin(), node_names_string.end(), back_inserter(node_names_string_sorted));
			this -> N_ = node_names_string_sorted.size();
		}
		assert(locked());
	}
};

NodeSet * network :: build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeNameType node_name_type) {
	ifstream edgelist(edgeListFileName);
	if(!edgelist) {
		cerr << "Error: edge list file (" << edgeListFileName << ") not found. Exiting" << endl;
		exit(1);
	}
	NodeSet * nodes = (node_name_type == NODE_NAME_INT64) ? static_cast<NodeSet*>(new NodeSet_<int64_t>) : new NodeSet_<string>;

	string line;
	int64_t line_num = 0;
	while (getline(edgelist, line)) {
		++ line_num;
		// PP(line);
		string left_node_string, right_node_string;
		istringstream fields(line);
		fields >> left_node_string;
		fields >> right_node_string;
		try {
			if(!fields)
				throw EdgeListFileFormatException();
			nodes->insert_string_version_of_name(left_node_string);
			nodes->insert_string_version_of_name(right_node_string);
		} catch (const EdgeListFileFormatException_ExpectedAnInt &) {
			cerr << "Error in edge list file on line " << line_num << ". Use the --str option if you wish to use non-integer node names. Exiting. : <" << line << ">" << endl;
			exit(1);
		} catch (const EdgeListFileFormatException &) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
	}
	nodes->finish_me();
	return nodes;
}
EdgeSet * network :: build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, NodeSet * node_set) {
	ifstream edgelist(edgeListFileName);
	if(!edgelist) {
		cerr << "Error: edge list file (" << edgeListFileName << ") not found. Exiting" << endl;
		exit(1);
	}

	// We don't need to think explicitly about int64_t here,
	// The file will give us strings, and we just need to map
	// those strings to a number between 0 and N.
	// In other words, the distinction between int64_t and string is
	// no longer really relevant.
	unordered_map<string, int> map_string_to_id;
	for(int i=0; i<node_set->N(); ++i) {
		map_string_to_id[ node_set->as_string(i) ] = i;
	}
	assert(node_set->N() == (int)map_string_to_id.size());

	EdgeSet * edges = new EdgeSet(weight_type);

	string line;
	int64_t line_num = 0;

	while (getline(edgelist, line)) {
		++ line_num;
		// PP(line);
		string left_node_string, right_node_string;
		istringstream fields(line);
		fields >> left_node_string;
		fields >> right_node_string;
		if(!fields) {
			cerr << "Error in edge list file on line " << line_num << ". Exiting. : <" << line << ">" << endl;
			exit(1);
		}
		int left_node_id, right_node_id;
		left_node_id = map_string_to_id[left_node_string];
		right_node_id = map_string_to_id[right_node_string];
		// If either of these two assertions fail, there is something weird going on.
		// I assumed the file was fine, due to the fact that build_node_set_from_edge_list
		// has already succeeded if we get in here
		assert(left_node_string == node_set->as_string(left_node_id));
		assert(right_node_string == node_set->as_string(right_node_id));
		EdgeSet :: Edge e;
		e.left = left_node_id;
		e.right = right_node_id;
		edges->edges.push_back(e);
	}
	assert(edges->edges.size() == line_num);
	return edges;
}
network :: Junction :: Junction(int edge_id_, int junction_type_, int this_node_id_, int far_node_id_)
			: edge_id(edge_id_), junction_type(junction_type_), this_node_id(this_node_id_), far_node_id(far_node_id_) {
}
bool network :: Junction :: operator <(const Junction &other) const {
	if(this->this_node_id < other.this_node_id) return true;
	if(this->this_node_id > other.this_node_id) return false;
	assert(this->this_node_id == other.this_node_id);
	if(this->far_node_id < other.far_node_id) return true;
	if(this->far_node_id > other.far_node_id) return false;
	assert(this->far_node_id == other.far_node_id);
	return false;
}
void network :: Junctions :: finish() {
	std :: sort(this->all_junctions_sorted.begin(), this->all_junctions_sorted.end());
}
network :: Junctions * network :: build_junctions_set_from_edges(const network :: EdgeSet * edge_set, const bool directed) {
	network :: Junctions * junctions = new network :: Junctions;
	const int E = edge_set->E();
	for(int e=0; e < E; ++e) {
		const network :: EdgeSet :: Edge & edge = edge_set->edges.at(e);
		network :: Junction l(e, directed ? 1 : 0, edge.left, edge.right);
		network :: Junction r(e, directed ?-1 : 0, edge.right, edge.left);
		junctions->all_junctions_sorted.push_back(l);
		junctions->all_junctions_sorted.push_back(r);
	}
	junctions->finish(); // tell junctions that we've finished telling it about the edges, it's ready to sort itself.
	assert(junctions->all_junctions_sorted.size() == (size_t)2*E);
	return junctions;
}
