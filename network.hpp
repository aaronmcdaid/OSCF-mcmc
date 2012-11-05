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
		int w_int;
	};
	std :: vector< Edge > edges;
};
EdgeSet * build_edge_set_from_edge_list(std :: string edgeListFileName, enum network :: EdgeSet :: WeightType weight_type, NodeSet_I * node_set);

} // namespace network
