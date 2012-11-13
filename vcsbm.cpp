const char gitstatus[] = 
#include "comment.txt"
#include "gitstatus.txt"
;
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
#include<algorithm>
#include<iomanip>
#include"format_flag_stack/format_flag_stack.hpp"


#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace std :: tr1;
using namespace network;

void vcsbm(const Network * net);

format_flag_stack :: FormatFlagStack stack;

int main(int argc, char **argv) {
	// Parse the args - there should be exactly one arg, the edge list
	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(1) ;
	if(args_info.git_version_flag) {
		PP(gitstatus);
		for (int i=0; i<argc; i++) {
			PP(argv[i]);
		}
	}
	if(args_info.inputs_num != 1) {
		cmdline_parser_print_help();
		exit(1);
	}

	const char * edgeListFileName   = args_info.inputs[0];

	const network :: Network * net = network :: build_network(edgeListFileName, args_info.stringIDs_flag, args_info.directed_flag, args_info.weighted_flag);

	vcsbm(net);
}

const int J = 10; // fix the upper bound on K at 10.
static void should_be_positive_(long double &x, const int line_no) {
	assert(x==x); // it's non NaN
	if(x<0.0L) {
		if(!VERYCLOSE(x,0.0L)) {
			PP2(line_no, x);
		}
		assert(VERYCLOSE(x,0.0L));
		x = 0.0L;
	}
}
#define SHOULD_BE_POSITIVE(x) should_be_positive_(x, __LINE__);


struct Q {
	// the variational lower bound is a function of Q.
	// We must keep track of various convenient summaries of Q, such as n_k,
	// but this will be managed in a separate class using callbacks.
	// This Q class will be kept very simple
	const int N;
	const int J;
	vector< vector<long double> > Q_;
	Q(const int N_, const int J_) : N(N_), J(J_), Q_(N, vector<long double>(J) ) {
	}
	long double get(int i, int k) const {
		return this->Q_.at(i).at(k);
	}
	struct Setter {
		Q * q;
		int i;
		int k;
		Setter(Q* q_, int i_, int k_) : q(q_), i(i_), k(k_) {}
		void operator = (long double val) const {
			this->q->set_one_cell(this->i, this->k, val);
		}
	};
	const Setter set (int i, int k) {
		return Setter(this, i, k);
	}
	void set_one_cell(int i, int k, long double val) {
		assert_0_to_1(val);
		try {
			const long double old_val = this->Q_.at(i).at(k);
			if(old_val == val)
				return;
			this->Q_.at(i).at(k) = val;
			this->notify_listeners(i,k,old_val, val);
		} catch (std :: out_of_range &e) {
			PP(this->Q_.size());
			PP(this->Q_.at(0).size());
			PP3(__LINE__, i, k);
			throw;
		}
	}

	// While the Q class itself is pretty simple,
	// there are other classes that will want to be notified
	// of any change to Q.  This is done via the Q_listener interface
	struct Q_listener {
		// any struct that wants to be notified of any changes to Q must
		// implement this interface.
		virtual void notify(int i, int k, long double old_val, long double new_val) = 0;
		virtual ~Q_listener() {}
	};
	vector<Q_listener *> listeners;
	void add_listener(Q_listener *l) {
		this->listeners.push_back(l);
	}
	void notify_listeners(int i,int k,long double old_val, long double new_val) {
		For(l, this->listeners) {
			(*l)->notify(i, k, old_val, new_val);
		}
	}
};

void dump(const Q *q, const Network *net) {
		cout << "node_id\t";
		cout << "node_name\t";
		cout << "degree\t";
		for(int k=0; k<J; ++k) {
			cout << "\tn_" << k;
		}
		cout << endl;
	for(int i=0; i<q->N; i++) {
		cout << i << '\t';
		cout << net->node_set->as_string(i) << '\t';
		cout << net->i.at(i).total_degree();
		for(int k=0; k<J; ++k) {
			const long double Q_ik = q->get(i,k);
			cout << '\t' << stack.push << fixed << setw(7) << setprecision(5) << Q_ik << stack.pop;
			// if(Q_ik > 0.5) cout << '+';
		}
		cout << endl;
	}
}

// Here are a few classes that will 'listen' to changes to Q
// and record various convenient summaries of Q, such as n_k
//
// mu_n_k = \sum_{i=1}^N Q_{ik}
// \sum_{i=1}^N Q_{ik}^2

template<int power>
struct Q_template_n_k : public Q :: Q_listener {
	vector<long double> n_k;
	Q_template_n_k() : n_k(10000) {
		assert(power == 1 || power == 2);
	}
	virtual void notify(int, int k, long double old_val, long double new_val) {
		if(power!=1) {
			assert(power==2);
			old_val = old_val * old_val;
			new_val = new_val * new_val;
		}
		assert_0_to_1(old_val);
		assert_0_to_1(new_val);
		this->n_k.at(k) -= old_val;
		if(this->n_k.at(k)<0.0L) {
			assert(VERYCLOSE(this->n_k.at(k) , 0.0L));
			this->n_k.at(k) = 0.0L;
		}
		this->n_k.at(k) += new_val;
	}
	void dump_me() const {
		cout << " n_k[0:20] : ";
		for(int k=0; k<20; ++k) {
			cout << ' ' << this->n_k.at(k);
		}
		cout << " ..." << endl;
	}
	virtual ~Q_template_n_k() {}
	void verify(const Q &q) const {
		vector<long double> verify_n_k(this->n_k.size());
		assert(!verify_n_k.empty());
		For(node, q.Q_) {
			for(int k=0; k < (int) node->size(); ++k) {
				long double Q_ik = node->at(k);
				if(power==1)
					verify_n_k.at(k) += Q_ik;
				else
					verify_n_k.at(k) += Q_ik * Q_ik;
			}
		}
		assert(this->n_k == verify_n_k);
	}
};
typedef Q_template_n_k<1> Q_mu_n_k;
typedef Q_template_n_k<2> Q_squared_n_k;

struct Q_entropy : public Q :: Q_listener {
	mutable long double entropy;
	Q_entropy() : entropy(0.0L) {}
	// this is          \EE ( - \log Q_ik )
	// this is  \int (- Q_ik \log Q_ik)
	// Note the minus, this is the entropy, not the negative entropy.
	virtual void notify(int, int, long double old_val, long double new_val) {
		assert_0_to_1(old_val);
		assert_0_to_1(new_val);
		if(old_val != 0.0L) {
			this->entropy -= - old_val * logl(old_val);
		}
		assert(this->entropy >= 0.0L);
		if(new_val != 0.0L) {
			this->entropy += - new_val * logl(new_val);
		}
		assert(this->entropy >= 0.0L);
	}
	void verify(const Q &q) const {
		long double verify_entropy = 0.0L;
		For(node, q.Q_) {
			For(cell, *node) {
				long double Q_ik = *cell;
				assert_0_to_1(Q_ik);
				if(Q_ik>0.0)
					verify_entropy += - Q_ik * logl(Q_ik);
			}
		}
		assert(VERYCLOSE(verify_entropy, this->entropy));
		this->entropy = verify_entropy;
	}
};

template<int power>
struct Q_templated_y_kl : public Q :: Q_listener {
	// Summing ONLY across the edges,
	// this will store Q_ik \times Q_jl,
	// or their square.
	// Must take care of direction too.
	// If directed, y_kl is the edges from k to l
	// If undirected, make k the smaller cluster-id
	unordered_map< pair<int,int> , long double> y_kl;
	const Network * const net;
	Q_templated_y_kl(const Network *net_) : net(net_) {}

	virtual void notify(int i, int , long double , long double ) {
		// Must consider all the Junctions at node i
		// .. but beware of self loops.
		For(junc_id, net->i.at(i).my_junctions) {
			const Junction & junc = net->junctions->all_junctions_sorted.at(*junc_id);
			assert(junc.this_node_id == i);
		}
	}
};
typedef Q_templated_y_kl<1> Q_mu_y_kl;
typedef Q_templated_y_kl<2> Q_squared_y_kl;

gsl_rng * global_r = NULL;

// A few globals (sorry) but they only help in verification and assertions
const Q_entropy *global_q_entropy;

// Hyperparameters
const long double alpha_for_stick_breaking = 1.0L;
const long double beta_1 = 1.0L;
const long double beta_2 = 1.0L;


static long double exp_log_Gamma_Normal(const long double mean, const long double variance) {
	assert(variance >= 0.0L);
	assert(mean >= 0.0L);
	return (mean - 0.5L) * logl(mean)
		+ (0.5L) * variance * ( 1.0L/mean + 0.25L/mean/mean )
		- mean
		+ 0.5L * logl(2 * M_PI);
}
void test_exp_log_Gamma_Normal(const long double mu, const long double stddev) {
	long double total = 0.0;
	const int S = 10000;
	for(int s = 0; s<S; ++s) {
		const long double random = mu + gsl_ran_gaussian(global_r, stddev);
		const long double x = log(gsl_sf_gamma(random));
		total += x;
	}
	PP4(mu, stddev, total/S, exp_log_Gamma_Normal(mu, stddev * stddev));
}
void test_exp_log_Gamma_Normal() {
	test_exp_log_Gamma_Normal(15,0.3);
	test_exp_log_Gamma_Normal(15,0.1);
	test_exp_log_Gamma_Normal(15,0.01);
	test_exp_log_Gamma_Normal(15,0.001);

	test_exp_log_Gamma_Normal(5,0.3);
	test_exp_log_Gamma_Normal(5,0.1);
	test_exp_log_Gamma_Normal(5,0.01);
	test_exp_log_Gamma_Normal(5,0.001);

	test_exp_log_Gamma_Normal(4,0.3);
	test_exp_log_Gamma_Normal(4,0.1);
	test_exp_log_Gamma_Normal(4,0.01);
	test_exp_log_Gamma_Normal(4,0.001);

	test_exp_log_Gamma_Normal(3,0.3);
	test_exp_log_Gamma_Normal(3,0.1);
	test_exp_log_Gamma_Normal(3,0.01);
	test_exp_log_Gamma_Normal(3,0.001);

	test_exp_log_Gamma_Normal(2,0.3);
	test_exp_log_Gamma_Normal(2,0.1);
	test_exp_log_Gamma_Normal(2,0.01);
	test_exp_log_Gamma_Normal(2,0.001);

	test_exp_log_Gamma_Normal(0.1,0.01);
	test_exp_log_Gamma_Normal(0.1,0.001);
	test_exp_log_Gamma_Normal(0.1,0.0001);
	test_exp_log_Gamma_Normal(0.1,0.00001);
}
static long double gamma_k(const int k) {
	assert(k>=0);
	assert(k<J);
	return powl(alpha_for_stick_breaking / (1.0L+alpha_for_stick_breaking), k);
}


long double calculate_first_four_terms_slowly(const Q *q, const Network * net, const bool everything_assigned_therefore_test = false) {
	const int N = q->N;
	const int E = net->edge_set->E();
	vector<long double> mu_n_k(J);
	vector<long double> sq_n_k(J);
	for(int i=0; i<N; ++i) {
		for(int k=0; k<J; ++k) {
			mu_n_k.at(k) += q->get(i,k);
			sq_n_k.at(k) += q->get(i,k)*q->get(i,k);
		}
	}

	if(everything_assigned_therefore_test) { // assert that the sum over mu_n_k == N
		long double  sum_of_mu_n_k = 0;
		For(x, mu_n_k) {
			sum_of_mu_n_k += *x;
		}
		assert(VERYCLOSE(sum_of_mu_n_k, N));
	}

	vector< vector<long double> > mu_y_kl(J, vector<long double>(J) );
	vector< vector<long double> > sq_y_kl(J, vector<long double>(J) );
	assert(E == (int)net->edge_set->edges.size());
	for(int m=0; m<E; ++m) {
		EdgeSet :: Edge edge = net->edge_set->edges.at(m);
		int i = edge.left;
		int j = edge.right;
		for(int k=0; k<J; ++k) {
			for(int l=0; l<J; ++l) {
				const long double Qik = q->get(i,k);
				const long double Qjl = q->get(j,l);
				int k2 = k;
				int l2 = l;
				{ // if undirected, we should have k2 <= j2
					if(k2>l2)
						swap(k2,l2);
				}
				mu_y_kl.at(k2).at(l2) += Qik * Qjl;
				sq_y_kl.at(k2).at(l2) += Qik * Qjl * Qik * Qjl;
			}
		}
	}

	// These next structures should be removed sometimes.  For now,
	// they are just for paranoid testing of p_kl
	vector< vector<long double> > mu_slowp_kl(J, vector<long double>(J) );
	vector< vector<long double> > sq_slowp_kl(J, vector<long double>(J) );
	for(int i=0; i<N; ++i) {
	for(int j=i; j<N; ++j) { // undirected j starts at i
		for(int k=0; k<J; ++k) {
			for(int l=0; l<J; ++l) {
				const long double Qik = q->get(i,k);
				const long double Qjl = q->get(j,l);
				int k2 = k;
				int l2 = l;
				{ // if undirected, we should have k2 <= j2
					if(k2>l2)
						swap(k2,l2);
				}
				mu_slowp_kl.at(k2).at(l2) += Qik * Qjl;
				sq_slowp_kl.at(k2).at(l2) += Qik * Qjl * Qik * Qjl;
			}
		}
	}
	}
	if(everything_assigned_therefore_test) {
		long double verify_num_pairs = 0.0L;
		for(int k=0; k<J; ++k) {
			for(int l=0; l<J; ++l) {
				if(l<k) {
					assert(0==mu_slowp_kl.at(k).at(l));
					assert(0==sq_slowp_kl.at(k).at(l));
					continue;
				}
				PP2(mu_slowp_kl.at(k).at(l), sq_slowp_kl.at(k).at(l));
				verify_num_pairs += mu_slowp_kl.at(k).at(l);
			}
		}
		assert(VERYCLOSE(verify_num_pairs , N * (N+1) / 2));
	}

	if(everything_assigned_therefore_test) { // assert that the sum over mu_y_kl == E
		long double sum_of_edge_partial_memberships = 0;
		For(one_cluster, mu_y_kl) {
			For(one_block, *one_cluster) {
				sum_of_edge_partial_memberships += *one_block;
			}
		}
		assert(VERYCLOSE(sum_of_edge_partial_memberships , E));
	}

	// Need to count the p_kl for the self loops.
	// .. maybe this can be ignored for the directed+selfloop model,
	// but otherwise it's important.
	vector< vector<long double> > mu_psl_kl(J, vector<long double>(J) );
	vector< vector<long double> > sq_psl_kl(J, vector<long double>(J) );
	for(int i=0; i<N; ++i) {
		for(int k=0; k<J; ++k) {
			for(int l=0; l<J; ++l) {
				const long double Qik = q->get(i,k);
				const long double Qil = q->get(i,l);
				mu_psl_kl.at(k).at(l) += Qik * Qil;
				sq_psl_kl.at(k).at(l) += Qik * Qil * Qik * Qil;
			}
		}
	}

	long double first_4_terms = 0.0L;

	// First term, E log Gamma( n_k + gamma_k )
	for(int k=0; k<J; k++) {
		const long double mu = mu_n_k.at(k);
		const long double var = mu - sq_n_k.at(k);
		assert(mu >= 0.0L);
		assert(var >= 0.0L);
		first_4_terms += exp_log_Gamma_Normal( mu + gamma_k(k), var );
	}
	// Second term, E log Gamma (y_kl + Beta_1)
//*
	// cout << endl;
	for(int k=0; k<J; k++) {
		for(int l=0; l<J; l++) {
			if(l<k) {
				// if undirected, these block don't really count.
				// In particular y_kl should be zero
				assert(0 == mu_y_kl.at(k).at(l));
				assert(0 == sq_y_kl.at(k).at(l));
				continue;
			}
			long double mu = mu_y_kl.at(k).at(l);
			long double var_y_kl = mu - sq_y_kl.at(k).at(l);

			long double mu_p_kl = mu_n_k.at(k)*mu_n_k.at(l);
			long double var_p_kl = mu_n_k.at(k)*mu_n_k.at(l) - sq_n_k.at(k) * sq_n_k.at(l);
			if(l==k) { // if undirected
				// Add in the self loops and then divide everything by two
				mu_p_kl  += mu_psl_kl.at(k).at(l);
				var_p_kl += mu_psl_kl.at(k).at(l) - sq_psl_kl.at(k).at(l);
				mu_p_kl /= 2;
				var_p_kl /= 2;
			} else { // bizarrelly, the self loops can contribute
			         // to between-cluster blocks.
				mu_p_kl  += mu_psl_kl.at(k).at(l);
				var_p_kl += mu_psl_kl.at(k).at(l) - sq_psl_kl.at(k).at(l);
			}
			const long double mu_slowp_KL = mu_slowp_kl.at(k).at(l);
			const long double var_slowp_KL = mu_slowp_kl.at(k).at(l) - sq_slowp_kl.at(k).at(l);

			//for the non-edges
			long double nonEdge_mu = mu_p_kl - mu;
			long double nonEdge_var = var_p_kl - var_y_kl;

			/*
			cout << k << ',' << l
				<< '\t' << mu << '(' << var_y_kl << ')'
				<< '\t' << mu_p_kl << '(' << var_p_kl << ')'
				<< '=' << mu_slowp_KL << '(' << var_slowp_KL << ')'
				<< '\t' << nonEdge_mu << '(' << nonEdge_var << ')'
				<< endl;
				*/

			SHOULD_BE_POSITIVE(mu);
			SHOULD_BE_POSITIVE(var_y_kl);
			SHOULD_BE_POSITIVE(mu_p_kl);
			SHOULD_BE_POSITIVE(var_p_kl);
			SHOULD_BE_POSITIVE(nonEdge_mu);
			SHOULD_BE_POSITIVE(nonEdge_var);
			assert(VERYCLOSE(mu_slowp_KL , mu_p_kl));
			assert(VERYCLOSE(var_slowp_KL , var_p_kl));

			first_4_terms += exp_log_Gamma_Normal( mu + beta_1, var_y_kl );
			first_4_terms += exp_log_Gamma_Normal( nonEdge_mu + beta_2, nonEdge_var );
			first_4_terms -= exp_log_Gamma_Normal( mu_p_kl + beta_1 + beta_2, var_p_kl );
		}
	}
// */

	// NOTE, we DO NOT include the entropy term in the return.

	return first_4_terms;
}
// The code above calculates stuff, below we have the actual algorithm.

void one_node_all_k(Q *q, const Network * net, const int node_id, const bool every_node_already_assigned = false /* *probably* not all assigned */);

void one_random_node_all_k(Q *q, const Network * net) {
	const int N = q->N;
	const double unif = gsl_rng_uniform(global_r);
	const int random_node = N * unif;
	one_node_all_k(q, net, random_node);
}

void one_node_all_k(Q *q, const Network * net, const int node_id, const bool every_node_already_assigned /* = false *probably* not all assigned */) {
	// pick one node at random,
	// remove it from all its clusters,
	// reassign to one cluster at a time
	const int num_clusters = q->Q_.at(node_id).size();
	assert(J==num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		q->set(node_id, k) = 0;
	}
	// store the baseline ?
	vector<long double> scores(J);
	for (int k = 0; k < num_clusters; ++k) {
		q->set(node_id, k) = 1;
		scores.at(k) = calculate_first_four_terms_slowly(q, net, every_node_already_assigned);
		q->set(node_id, k) = 0;
	}
	const long double max_score = *max_element(scores.begin(), scores.end());
	For(score, scores) {
		// PP(*score);
		*score -= max_score;
	}
	// For(score, scores) { PP(*score); }
	long double total = 0.0L;
	For(score, scores) {
		total += expl(*score);
	}
	For(score, scores) {
		*score  = expl(*score) / total;
	}
	long double check_new_total_is_1 = 0.0L;
	For(score, scores) {
		// PP(*score);
		check_new_total_is_1 += *score;
	}
	// PP(check_new_total_is_1 - 1.0L);
	assert(VERYCLOSE(check_new_total_is_1, 1.0L));
	// scores.front() += 1.0L - check_new_total_is_1;
	// assert(scores.front()>=0.0L);

	for(int k=0; k<J; ++k) {
		q->set(node_id, k) = scores.at(k);
	}
}

void vcsbm(const Network * net) {
	const int N = net->N();

	global_r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(global_r, args_info.seed_arg);

	// test_exp_log_Gamma_Normal(); return;

	Q q(N,J);
	Q_mu_n_k mu_n_k;
	Q_squared_n_k squared_n_k;
	Q_entropy entropy;
	Q_squared_y_kl squared_y_kl(net);

	PP(entropy.entropy);

	q.add_listener(&mu_n_k);
	q.add_listener(&squared_n_k);
	q.add_listener(&entropy);
	q.add_listener(&squared_y_kl);
	global_q_entropy = &entropy;


	entropy.verify(q);
	mu_n_k.dump_me();
	squared_n_k.dump_me();
	PP(entropy.entropy);
	entropy.verify(q);
	mu_n_k.verify(q);
	squared_n_k.verify(q);

	for(int i=0; i<N; i++) {
		one_node_all_k(&q, net, i);
		PP(entropy.entropy);
		mu_n_k.dump_me();
	}
	// everything assigned somewhere
	cout << "everything assigned somewhere" << endl;
	// calculate_first_four_terms_slowly(&q, net, true);
	for(int repeat = 0; repeat < 3; ++repeat) {
		PP(repeat);
		for(int i=0; i<N; i++) {
			one_node_all_k(&q, net, i, true);
			PP(entropy.entropy);
			mu_n_k.dump_me();
		}
	}
	dump(&q);
}
