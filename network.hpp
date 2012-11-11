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
	virtual std :: string as_string(int) const = 0;
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

const NodeSet * build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeNameType node_name_type);

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
	inline int E() const { return this->edges.size(); }
};
const EdgeSet * build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, const NodeSet * node_set);

struct Junction {
	int edge_id; // edge_id: a number between 0 and E
	int junction_type; // 1: I'm the source, and the other end of the sink
	                   // -1: Other way around
			   // 0: Undirected, so it doesn't matter.
	int this_node_id; // id of this node at this end of the edge.
	int far_node_id; // id of the node at the other end of this  edge.
	Junction(int edge_id_, int junction_type_, int this_node_id_, int far_node_id_);
	bool operator <(const Junction &other) const;
};
struct Junctions {
	std :: vector<Junction> all_junctions_sorted; // we'll sort this later
	void finish();
};
const Junctions * build_junctions_set_from_edges(const EdgeSet * edge_set, bool directed);

} // namespace network
