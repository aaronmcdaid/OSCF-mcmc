#include<vector>
namespace network {

struct NodeSet_I { // the interface to a set of nodes. Either int64, or strings
	enum NodeNameType { NODE_NAME_INT64, NODE_NAME_STRING };
	virtual NodeNameType get_NodeNameType() const = 0;
	virtual int N() const = 0;
	virtual std :: string as_string(int node_id) const = 0;
	virtual int64_t as_int64(int node_id) const = 0;
	virtual ~NodeSet_I() {
	}
};

NodeSet_I * build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeSet_I :: NodeNameType node_name_type);

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
EdgeSet * build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, NodeSet_I * node_set);

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
