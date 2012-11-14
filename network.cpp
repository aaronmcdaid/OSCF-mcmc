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
	virtual string as_string(int) const;
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
	string NodeSet_<string> :: as_string(int node_id) const {
		assert(node_id >= 0);
		assert(node_id < (int)this->vector_of_names.size());
		return this->vector_of_names.at(node_id);
	}
	template<>
	string NodeSet_<int64_t> :: as_string(int node_id) const {
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

const NodeSet * network :: build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeNameType node_name_type) {
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
static const EdgeSet * build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, const NodeSet * node_set, const bool directed_flag) {
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
	{
		// I should warn about duplicate edges
		set< pair<int,int> > unique_edges;
		bool seen_duplicates_already = false;
		For(edge, edges->edges) {
			int n1 = edge->left;
			int n2 = edge->right;
			if(!directed_flag) {
				if(n1 > n2) {
					swap(n1,n2);
				}
				assert(n1<=n2);
			}
			const bool was_inserted = unique_edges.insert( make_pair(n1,n2) ).second;
			if(!was_inserted && !seen_duplicates_already) {
				cout << "Duplicate edge(s). Discarding it as I don't support weights yet." << endl;
				seen_duplicates_already = true;
			}
		}
		unless(unique_edges.size() == edges->edges.size()) {
			edges->edges.clear();
			For(node_pair, unique_edges) {
				EdgeSet :: Edge e;
				e.left = node_pair->first;
				e.right = node_pair->second;
				edges->edges.push_back(e);
			}
		}
		assert(unique_edges.size() == edges->edges.size());
	}
	return edges;
}
network :: Junction :: Junction(int edge_id_, int junction_type_, int this_node_id_, int far_node_id_, bool i_am_the_second_self_loop_junction_)
			: edge_id(edge_id_), junction_type(junction_type_), this_node_id(this_node_id_), far_node_id(far_node_id_), i_am_the_second_self_loop_junction(i_am_the_second_self_loop_junction_) {
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
const network :: Junctions * network :: build_junctions_set_from_edges(const network :: EdgeSet * edge_set, const bool directed) {
	network :: Junctions * junctions = new network :: Junctions;
	const int E = edge_set->E();
	for(int e=0; e < E; ++e) {
		const network :: EdgeSet :: Edge & edge = edge_set->edges.at(e);
		network :: Junction l(e, directed ? 1 : 0, edge.left, edge.right, false);
		network :: Junction r(e, directed ?-1 : 0, edge.right, edge.left, edge.right == edge.left);
		junctions->all_junctions_sorted.push_back(l);
		junctions->all_junctions_sorted.push_back(r);
	}
	junctions->finish(); // tell junctions that we've finished telling it about the edges, it's ready to sort itself.
	assert(junctions->all_junctions_sorted.size() == (size_t)2*E);
	return junctions;
}
const network :: Network * network :: build_network(const std :: string file_name, const bool stringIDs_flag, const bool directed_flag, const bool weighted_flag) {
	assert(weighted_flag == false);
	const network :: NodeSet * node_set = build_node_set_from_edge_list(file_name,
			stringIDs_flag
			?  network :: NODE_NAME_STRING
			:  network :: NODE_NAME_INT64
			);
	const network :: EdgeSet * edge_set = build_edge_set_from_edge_list(file_name,
			weighted_flag
			? network :: EdgeSet :: WEIGHT_INT
			: network :: EdgeSet :: WEIGHT_NONE
			, node_set
			, directed_flag
			);
	const int N = node_set -> N();
	const int E = edge_set->edges.size();
	PP2(N,E);

	// Finally, we create the list of Junctions - by doubling the list of Edges
	const network :: Junctions * junctions = build_junctions_set_from_edges(edge_set, directed_flag);
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
	Network * net = new Network;
	// I should be these next few lines into the Network constructor
	net->node_set = node_set;
	net->edge_set = edge_set;
	net->junctions = junctions;
	assert(net->i.empty());
	net->i.resize(net->N());
	for(int junc = 0; junc < (int)net->junctions->all_junctions_sorted.size(); ++junc) {
		const Junction & jun = net->junctions->all_junctions_sorted.at(junc);
		net->i.at(jun.this_node_id).my_junctions.push_back(junc);
	}
	return net;
}
