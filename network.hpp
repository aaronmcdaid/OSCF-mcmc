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

struct Junction {
	int edge_id; // edge_id: a number between 0 and E
	int junction_type; // 1: I'm the source, and the other end of the sink
	                   // -1: Other way around
			   // 0: Undirected, so it doesn't matter.
	int this_node_id; // id of this node at this end of the edge.
	int far_node_id; // id of the node at the other end of this  edge.
	bool i_am_the_second_self_loop_junction;
	Junction(int edge_id_, int junction_type_, int this_node_id_, int far_node_id_, bool i_am_the_second_self_loop_junction);
	bool operator <(const Junction &other) const;
};
struct Junctions {
	std :: vector<Junction> all_junctions_sorted; // we'll sort this later
	void finish();
};
const Junctions * build_junctions_set_from_edges(const EdgeSet * edge_set, bool directed);

struct OneNode {
	std :: vector<int> my_junctions; // all of the junction ids associated with this node
	inline int total_degree() const { return this->my_junctions.size(); }
};

struct Network__unConst { // bring all the above together
	const network :: NodeSet * node_set;
	const network :: EdgeSet * edge_set;
	const network :: Junctions * junctions;
	std :: vector<OneNode> i;
	inline int N() const { return this->node_set->N(); }
	inline int E() const { return this->edge_set->E(); }
};
typedef const Network__unConst Network; // The exposed Network should always be const.


const Network * build_network(std :: string file_name, const bool stringIDs_flag, const bool directed_flag, const bool weighted_flag);

} // namespace network
