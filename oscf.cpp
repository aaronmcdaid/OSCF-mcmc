#include "gitstatus.hpp"
#include "cmdline.h"
gengetopt_args_info args_info; // a global variable! Sorry.
#include "macros.hpp"

#include "network.hpp"

#include<cassert>
#include<stdexcept>
#include<cmath>
#include<tr1/unordered_map>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_sf_gamma.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_statistics_double.h>
#include<algorithm>
#include<iomanip>
#include<limits>
#include<algorithm>
#include<fstream>
#include<array>

//#include"format_flag_stack/format_flag_stack.hpp"

#include"state.hpp"
#include"score.hpp"
#include"moves.hpp"


#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace std :: tr1;
using namespace network;


// format_flag_stack :: FormatFlagStack stack;


void oscf(Net net);

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

	const network :: Network * net = network :: build_network(edgeListFileName, args_info.stringIDs_flag, args_info.directed_flag, args_info.weighted_flag);

	// if(args_info.GT_vector_arg) { loadGroundTruth(args_info.GT_vector_arg); assert(global_groundTruth.size() == net->N()); }

	//
	PP2(net->N(), net->E());
	
	seed_the_random_number_generator(args_info.seed_arg);

	oscf(net);
}

void dump_all(const State & st) {
	cout << endl << " ===" << endl;
	for(int k=0; k<st.get_K(); ++k) {
		cout << '\t';
		st.get_comms().at(k).dump_me();
	}
	cout << endl;
	cout << " ==" << endl << endl;
}
#define CHECK_PMF_TRACKER(track, actual) do { const long double _actual = (actual); long double & _track = (track); if(VERYCLOSE(_track,_actual)) { track = _actual; } else { PP(_actual - track); } assert(_track == _actual); } while(0)
void oscf(Net net) {
	State st(net); // initialize with every edge in its own community
	Score sc(st);
	long double cmf_track = sc.score();
	assert(st.every_edge_non_empty());
	PP(cmf_track);
	CHECK_PMF_TRACKER(cmf_track, sc.score());
	for (int rep = 0; rep < 250000; ++rep) {
		if(rep % 10000 == 0)
			cerr << rep << endl;
		if(rep % 1000 == 0)
			CHECK_PMF_TRACKER(cmf_track, sc.score());
		for(int64_t e = 0; e<net->E(); ++e) {
			cmf_track += gibbsUpdate(e, sc);
			cmf_track += metroK(sc);
		}
		dump_all(st);
	}
	assert(st.every_edge_non_empty());
	PP(cmf_track);
	CHECK_PMF_TRACKER(cmf_track, sc.score());
}


