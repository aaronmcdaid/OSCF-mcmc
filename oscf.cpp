#include "gitstatus.hpp"

#include<cassert>
#include<stdexcept>
#include<cmath>
#include<tr1/unordered_map>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_statistics_double.h>
#include<iomanip>
#include<limits>
#include<algorithm>
#include<fstream>
#include<sstream>

#include"lvalue_input.hpp"
#include"state.hpp"
#include"score.hpp"
#include"moves.hpp"
#include"onmi.hpp"
#include "cmdline.h"
gengetopt_args_info args_info; // a global variable! Sorry.
#include "macros.hpp"
#include "network.hpp"
#include "format_flag_stack/format_flag_stack.hpp"
static format_flag_stack :: FormatFlagStack stack;


#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace std :: tr1;
using namespace network;
using namespace lvalue_input;


// format_flag_stack :: FormatFlagStack stack;


void oscf(Net net, const gengetopt_args_info &args_info);

int main(int argc, char **argv) {
	// Parse the args - there should be exactly one arg, the edge list
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(1) ;
	if(args_info.git_version_flag) {
		cout << gitstatus;
		for (int i=0; i<argc; i++) {
			PP(argv[i]);
		}
		cout << "=====" << endl;
	}
	if(args_info.inputs_num != 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	const char * edgeListFileName   = args_info.inputs[0];

	const network :: Network * net = network :: build_network(edgeListFileName, args_info.stringIDs_flag, false/*args_info.directed_flag*/, false/*args_info.weighted_flag*/);

	// if(args_info.GT_vector_arg) { loadGroundTruth(args_info.GT_vector_arg); assert(global_groundTruth.size() == net->N()); }

	//
	PP2(net->N(), net->E());
	
	seed_the_random_number_generator(args_info.seed_arg);

	oscf(net, args_info);
}

static long double entropy_of_this_state(in< State > st) {
	long double entropy = 0.0;
	For(comm, st->get_comms()) {
		const int64_t num_u = comm->get_num_unique_nodes_in_this_community();
		if(num_u > 0) {
			const long double p_k = (long double) num_u / (long double) st->N;
			entropy += -p_k * log2l(p_k);
		}
	}
	return entropy;
}
static	bool	biggest_cluster_first(const OneCommunitySummary &l, const OneCommunitySummary &r) {
			if(l.num_unique_nodes_in_this_community == r.num_unique_nodes_in_this_community)
				return l.num_edges > r.num_edges;
			return l.num_unique_nodes_in_this_community > r.num_unique_nodes_in_this_community;
}
void dump_all(const State & st, const int64_t rep, in< std::vector< std::vector<int64_t> > > ground_truth, const bool /*cluster_sizes*/ /*= false*/) {
	cout << " ===" << endl;
	const long double entropy = entropy_of_this_state(st);
	const long double onmi = ground_truth->empty() ? -1.0L : calculate_oNMI(ground_truth, st);
	PP2(rep, ELAPSED());
	const int64_t true_K = ground_truth->empty() ? -1 : ground_truth->size();
	PP5(rep, entropy, st.get_K(), true_K, onmi);
	{
		cout << "average assignments per edge: " << double(st.total_count_of_edge_assignments) / st.E << '\t';
		int64_t printed_so_far = 0;
		assert(0 == st.frequencies_of_edge_occupancy.at(0));
		for(size_t f = 1; f<st.frequencies_of_edge_occupancy.size(); ++f) {
			ostringstream oss;
			const int64_t freq_f = st.frequencies_of_edge_occupancy.at(f);
			if(freq_f != 0)
				oss << st.frequencies_of_edge_occupancy.at(f) << '(' << f << ')';
			cout << setw(8) << oss.str();
			printed_so_far += f * freq_f;
			if(printed_so_far >= st.total_count_of_edge_assignments) {
				assert(printed_so_far == st.total_count_of_edge_assignments);
				break;
				cout << '*';
			}
		}
		cout << endl;
	}
	//if(cluster_sizes)
	{ // print sizes of all clusters?
		cout << st.get_K();
		vector< OneCommunitySummary > summaries; // we'll reorder them by size later
		for(int k=0; k<st.get_K(); ++k) {
			summaries.push_back( st.get_comms().at(k).get_one_community_summary());
		}
		const size_t how_many_to_print = 15;
		if (summaries.size() >= how_many_to_print)
			partial_sort(summaries.begin(), summaries.begin()+how_many_to_print, summaries.end(), biggest_cluster_first );
		else
			sort(summaries.begin(), summaries.end(), biggest_cluster_first );
		for(size_t k=0; (int)k<st.get_K(); ++k) {
			if(k>=how_many_to_print) {
				cout << "  ... too many ...";
				break;
			}
			cout << ";  ";
			in< OneCommunitySummary > comm_k = summaries.at(k);
			cout
				<< "  " << comm_k->num_unique_nodes_in_this_community
				<< stack.push << fixed <<  setprecision(2)
				<< ' '
				<< comm_k->density()
				<< stack.pop;
				;
		}
		cout << endl;
	}
	cout << " ==" << endl << endl;
}
void dump_truncated_node_cover(const State & st) {
	if(! args_info.save_z_arg || ! args_info.save_z_arg[0])
		return;
	ofstream z_file(args_info.save_z_arg);
	assert(args_info.save_z_arg && args_info.save_z_arg[0]);
	for(int k=0; k<st.get_K(); ++k) {
		vector<int64_t> nodes_in_this_comm = st.get_comms().at(k).get_my_nodes_NO_COUNT();
		sort(nodes_in_this_comm.begin(), nodes_in_this_comm.end());
		For(n, nodes_in_this_comm) {
			if(n!=nodes_in_this_comm.begin())
				z_file << ' ';
			string node_name = st.net->node_set->as_string(*n);
			z_file << node_name;
		}
		z_file << endl;
	}
}

static vector< vector<int64_t> > load_ground_truth(const NodeSet * const node_set, const char * gt_file_name) {
	unordered_map<string, int> map_string_to_id;
	for(int i=0; i<node_set->N(); ++i) {
		map_string_to_id[ node_set->as_string(i) ] = i;
	}
	assert(node_set->N() == (int)map_string_to_id.size());

	vector< vector<int64_t> > all_gt_comms;

	ifstream gt(gt_file_name);
	string line;
	while(getline(gt, line)) {
		vector<int64_t> this_comm;
		istringstream fields(line);
		string field;
		while(fields >> field) {
			const int node_id = map_string_to_id[field];
			this_comm.push_back(node_id);
			assert(node_set->N() == (int)map_string_to_id.size());
		}
		assert(!this_comm.empty());
		all_gt_comms.push_back(this_comm);

	}

	// Finally, ensure the map hasn't increased in size!  This would mean that the GT has invalid entries in it.
	assert(node_set->N() == (int)map_string_to_id.size());
	assert(!all_gt_comms.empty());
	return all_gt_comms;
}

#define CHECK_PMF_TRACKER(track, actual) do { const long double _actual = (actual); long double & _track = (track); if(VERYCLOSE(_track,_actual)) { track = _actual; } else { PP(_actual - track); } assert(_track == _actual); } while(0)
void oscf(Net net, const gengetopt_args_info &args_info) {
	State st(net);
	Score sc(st);
	sc.m_iidBernoulli_arg = args_info.m_iidBernoulli_arg;
	assert(st.get_K() == 0);
	if(0) { // every edge in its own cluster
		for(int64_t e = 0; e < st.E; ++e) {
			const int64_t new_cluster_id = st.append_empty_cluster();
			assert(new_cluster_id == e);
			sc.add_edge(e, e);
		}
	}
	for(int identical_clusters=0; identical_clusters<1; ++identical_clusters) {
		const int64_t new_cluster_id = st.append_empty_cluster();
		//assert(0 == new_cluster_id);
		for(int64_t e = 0; e < st.E; ++e) {
			st.add_edge(e, new_cluster_id);
		}
	}
	long double cmf_track = sc.score();
	PP(cmf_track);
	if(args_info.init_seedExpand_flag) { // initialize with many runs of expand_seed
		for(int i = 0; i < 10000; ++i) {
			const int k = 0; // always expand from the initial "big" cluster
			if(sc.state.get_one_community_summary(k).num_edges == 0)
				break;
			const int qplus1 = sc.state.get_K();
			cmf_track += sc.append_empty_cluster();
			assert(k!=qplus1);
			vector<int> E = my_edges_in_a_random_order(sc.state, k); // I don't really need them to be random, but what the hell
			PP(i, E.size());
			const size_t seed_edge = E.front();
			expand_seed(k, qplus1, seed_edge, E, sc.state.net, cmf_track, sc);
			PP(sc.state.get_one_community_summary(qplus1).num_edges);
			if(
				sc.state.get_one_community_summary(qplus1).num_edges > sc.state.get_one_community_summary(0).num_edges
				||
				sc.state.get_one_community_summary(qplus1).num_edges == 0
			)
				break;
		}
	}
	PP(st.get_K());
	const bool K_can_vary = args_info.K_arg == -1;
	if( !K_can_vary ) {
		while(st.get_K() < args_info.K_arg) {
			cmf_track += sc.append_empty_cluster();
		}
		assert(st.get_K() == args_info.K_arg);
	}

	vector< vector<int64_t> > ground_truth; // leave empty if no ground truth was specified
	if(args_info.GT_arg) { //Load a ground Truth?
		ground_truth = load_ground_truth(net->node_set, args_info.GT_arg);
	}


	// Next, we ensure that each edge is fully assigned.
	// Do the test twice,
	assert(0 == st.get_frequencies_of_edge_occupancy().at(0));
	assert(st.every_edge_non_empty());

	PP(cmf_track);
	CHECK_PMF_TRACKER(cmf_track, sc.score());

	vector<int64_t> edges_in_random_order; // will randomize later, at the start of each iteration
	for(int64_t e = 0; e<net->E(); ++e) { edges_in_random_order.push_back(e); }
	assert((int64_t)edges_in_random_order.size() == net->E());

	vector<int64_t> nodes_in_random_order; // will randomize later, at the start of each iteration
	for(int64_t n = 0; n<net->N(); ++n) { nodes_in_random_order.push_back(n); }
	assert((int64_t)nodes_in_random_order.size() == net->N());

	//long double min_K_double = net->N() < net->E() ? net->N() : net->E();
	long double min_K_double = net->E();

	for (int rep = 1; rep <= args_info.iterations_arg; ++rep) {
		{ // check for the -K arg
			if(!K_can_vary)
				assert(st.get_K() == args_info.K_arg);
		}
		{
			min_K_double *= 0.97;
			if(min_K_double < 1 || rep >= 200)
				min_K_double = 1;
		}
		GLOBAL_constraint_min_K = 1; // min_K_double;
		if(GLOBAL_constraint_min_K>1)
			PP(GLOBAL_constraint_min_K);
		assert(GLOBAL_constraint_min_K >= 1);
		assert(GLOBAL_constraint_min_K <= st.get_K());

		random_shuffle(edges_in_random_order.begin(), edges_in_random_order.end());
		random_shuffle(nodes_in_random_order.begin(), nodes_in_random_order.end());
		size_t node_offset = 0;
		//if(!args_info.algo_seedSplit_flag)
		For(e, edges_in_random_order) {
			assert(st.get_K() >=  GLOBAL_constraint_min_K);
			if(args_info.metroK_algo_arg && K_can_vary) cmf_track += metroK(sc);
			assert(st.get_K() >=  GLOBAL_constraint_min_K);
			if(args_info.metro1Comm1Edge_algo_arg) cmf_track += gibbs_one_comm_one_edge(sc, *e);
			//cmf_track += gibbsUpdate(*e, sc);
			//cmf_track += one_node_simple_update(sc);
			if(args_info.NearbyGibbs_algo_arg) cmf_track += gibbsUpdateNearby(sc, *e);
			{ // every time we do an edge, do a node aswell
				if(args_info.Simplest1Node_algo_arg) cmf_track += one_node_SIMPLEST_update(sc, nodes_in_random_order.at(node_offset) );
				//PP2(node_offset, nodes_in_random_order.at(node_offset) );
				++node_offset;
				if(node_offset >= nodes_in_random_order.size())
					node_offset = 0;
			}

			if(gsl_ran_bernoulli(rng(), 0.05) && args_info.AnySM_algo_arg    && K_can_vary) cmf_track += split_or_merge(sc);
			if(gsl_ran_bernoulli(rng(), 0.05) && args_info.SharedSM_algo_arg && K_can_vary) cmf_track += split_or_merge_on_a_shared_edge(sc);
			if(gsl_ran_bernoulli(rng(), 0.05) && args_info.M3_algo_arg)                     cmf_track += M3(sc);
		}
		if(  //gsl_ran_bernoulli(rng(), 0.05) &&
				args_info.algo_seedSplit_flag && K_can_vary) cmf_track += split_or_merge_by_seed_expansion(sc);
			/*
		if(args_info.algo_seedSplit_flag) {
			//for(size_t i=0; i<10000; ++i) {
				cmf_track += split_or_merge_by_seed_expansion(sc);
				//CHECK_PMF_TRACKER(cmf_track, sc.score());
			//}
			//cout << "exiting early" << endl;
			//exit(1);
			For(e, edges_in_random_order) {
				assert(st.get_K() >=  GLOBAL_constraint_min_K);
				if(args_info.metroK_algo_arg && K_can_vary) cmf_track += metroK(sc);
				if(args_info.NearbyGibbs_algo_arg) cmf_track += gibbsUpdateNearby(sc, *e);
			}
		}
		*/
		/*
		for(int i=0; i<net->E(); ++i) {
			// if(i % 1000 == 0) cerr << i << ',' << st.get_K() << endl;
			if(args_info.AnySM_algo_arg    && K_can_vary) cmf_track += split_or_merge(sc);
			if(args_info.SharedSM_algo_arg && K_can_vary) cmf_track += split_or_merge_on_a_shared_edge(sc);
			if(args_info.M3_algo_arg)                     cmf_track += M3(sc);
		}
		*/

		if(rep>0 && rep % 10 == 0) {
			//cerr << rep << endl;
			assert(st.every_edge_non_empty());
			CHECK_PMF_TRACKER(cmf_track, sc.score());
		}
		if(rep>0 && rep % 1 == 0) { // || rep < 100) {
			dump_all(st, rep, ground_truth);
			dump_truncated_node_cover(st);
		}
	}
	PP2("fin", ELAPSED());
	assert(st.every_edge_non_empty());
	PP(cmf_track);
	CHECK_PMF_TRACKER(cmf_track, sc.score());
}


