#include<vector>
namespace network {

enum NodeNameType { NODE_NAME_INT64, NODE_NAME_STRING };
struct NodeSet {
	// A Nodeset will take a load of strings and store them (or their int equivalent)
	// Then, it will sort them (either lexicographically or numerically) and allow
	// you random access to the sorted vector
	//
	// So yes, the interface is entirely string-based, as this is all that you will need
	// most of the time.  But under the hood they may be stored, and sorted, as strings.
	virtual void insert_string_version_of_name(std :: string) = 0;
	virtual std :: string as_string(int) = 0;
	virtual void finish_me() = 0;
	virtual int N() const = 0;
	virtual ~NodeSet() {}
	virtual NodeNameType get_NodeNameType() const = 0;
};

struct NodeSet_I { // the interface to a set of nodes. Either int64, or strings
	virtual NodeNameType get_NodeNameType() const = 0;
	virtual int N() const = 0;
	virtual std :: string as_string(int node_id) const = 0;
	virtual int64_t as_int64(int node_id) const = 0;
	virtual ~NodeSet_I() {
	}
};

NodeSet * build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeNameType node_name_type);

struct EdgeSet {
	// This is a pretty dumb object, it just knows the edges in the order in which they appeared in the text file.
	// It does understand weights though.
	enum WeightType { WEIGHT_NONE, WEIGHT_INT };
	const WeightType weight_type;
	EdgeSet(WeightType weight_type_) : weight_type(weight_type_) {
	}
	struct Edge {
		int left;
		int right;
		WeightType w_type;
		union {
			int w_int;
		} w_union;
	};
	std :: vector< Edge > edges;
};
EdgeSet * build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, NodeSet * node_set);

struct EndPoint {
	int edge_id; // edge_id: a number between 0 and E
	int endpoint_type; // 1: I'm the source, and the other end of the sink
	                   // -1: Other way around
			   // 0: Undirected, so it doesn't matter.
	int this_node_id; // id of this node at this end of the edge.
	int far_node_id; // id of the node at the other end of this  edge.
	EndPoint(int edge_id_, int endpoint_type_, int this_node_id_, int far_node_id_);
	bool operator <(const EndPoint &other) const;
};
struct EndPoints {
	std :: vector<EndPoint> all_endpoints_sorted; // we'll sort this later
	void finish();
};

} // namespace network
