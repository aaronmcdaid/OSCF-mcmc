namespace network {

struct NodeSet_I { // the interface to a set of nodes. Either int64, or strings
	enum NodeNameType { NODE_NAME_INT64, NODE_NAME_STRING };
	virtual int get_NodeNameType() const = 0;
	virtual int N() const = 0;
	virtual std :: string as_string(int node_id) const = 0;
	virtual ~NodeSet_I() {
	}
};

NodeSet_I * build_node_set_from_edge_list(std :: string edgeListFileName, enum network :: NodeSet_I :: NodeNameType node_name_type);

} // namespace network
