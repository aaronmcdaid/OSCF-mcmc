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

#include"format_flag_stack/format_flag_stack.hpp"


#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace std :: tr1;
using namespace network;

void vcsbm(Network * net);

format_flag_stack :: FormatFlagStack stack;

vector<int> global_groundTruth;
vector<int> *global_groundTruth_ptr = NULL;
void loadGroundTruth(const char *filename) {
	assert(filename);
	ifstream file(filename);
	int line_count = 0;
	int z_i;
	while(file >> z_i) {
		++line_count;
		global_groundTruth.push_back(z_i);
	}
	global_groundTruth_ptr = &global_groundTruth;
}

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

	if(args_info.GT_vector_arg) {
		loadGroundTruth(args_info.GT_vector_arg);
		assert(global_groundTruth.size() == net->N());
	}

	vcsbm(net);
}

typedef vector< vector<long double> > VVL;

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

#define FIXED_K 7
const int J = FIXED_K; // fix the upper bound on K at 10.


struct Q {
	// the variational lower bound is a function of Q.
	// We must keep track of various convenient summaries of Q, such as n_k,
	// but this will be managed in a separate class using callbacks.
	// This Q class will be kept very simple
	const int N;
	const int J;
	VVL Q_;
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
		virtual void notify(int i, int k, long double old_val, long double new_val, const Q * q) = 0;
		virtual ~Q_listener() {}
	};
	vector<Q_listener *> listeners;
	void add_listener(Q_listener *l) {
		this->listeners.push_back(l);
	}
	void notify_listeners(int i,int k,long double old_val, long double new_val) {
		For(l, this->listeners) {
			(*l)->notify(i, k, old_val, new_val, this);
		}
	}
};

/*
 * We will have a number of Q_listener objects, tracking
 * useful summaries of the data in Q.
 *
 * There will also be a single Tracker object, which
 * will store all the Q_listener objects.
 * We can ask the Tracker for the score every now
 * and then, and also ask it to force a verification
 * of *everything*
 */

void dump(const Q *q, Network *net) {
		cout << "node_id\t";
		cout << "nodename";
		cout << " degree\t";
		for(int k=0; k<J; ++k) {
			cout << "     n_" << k;
		}
		cout << endl;
	for(int i=0; i<q->N; i++) {
		cout << i << '\t';
		cout << net->node_set->as_string(i) << '\t';
		cout << net->i.at(i).total_degree() << '\t';
		for(int k=0; k<J; ++k) {
			const long double Q_ik = q->get(i,k);
			cout << stack.push << fixed << setw(8) << setprecision(5) << Q_ik << stack.pop;
			// if(Q_ik > 0.5) cout << '+';
		}
		cout << endl;
	}
}

template<int Power>
static inline long double power(long double x);
template<>
inline long double power<1>(long double x) { return x; }
template<>
inline long double power<2>(long double x) { return x*x; }

// Here are a few classes that will 'listen' to changes to Q
// and record various convenient summaries of Q, such as n_k
//
// mu_n_k = \sum_{i=1}^N Q_{ik}
// \sum_{i=1}^N Q_{ik}^2

template<int Power>
struct Q_template_n_k : public Q :: Q_listener {
	mutable vector<long double> n_k;
	Q_template_n_k() : n_k(10000) { }
	virtual void notify(int, int k, long double old_val, long double new_val, const Q *) {
		old_val = power<Power>(old_val);
		new_val = power<Power>(new_val);
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
		cout << " n_k[0:J]\t\t";
		for(int k=0; k<J; ++k) {
			cout << stack.push << fixed << setw(8) << setprecision(3)
				<< this->n_k.at(k)
				<< stack.pop;
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
				verify_n_k.at(k) += power<Power>(Q_ik);
			}
		}
		assert(this->n_k.size() == verify_n_k.size());
		for(size_t k=0; k<this->n_k.size(); ++k) {
			assert(VERYCLOSE( this->n_k.at(k),   verify_n_k.at(k) ));
		}
		this->n_k = verify_n_k;
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
	virtual void notify(int, int, long double old_val, long double new_val, const Q *) {
		assert_0_to_1(old_val);
		assert_0_to_1(new_val);
		if(old_val != 0.0L) {
			this->entropy -= - old_val * logl(old_val);
		}
		SHOULD_BE_POSITIVE(this->entropy);
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
struct Q_sum_of_mu_n_k : public Q :: Q_listener {
	mutable long double sum_of_mu_n_k;
	Q_sum_of_mu_n_k() : sum_of_mu_n_k(0.0L) {}
	virtual void notify(int, int, long double old_val, long double new_val, const Q *) {
		this->sum_of_mu_n_k -= old_val;
		SHOULD_BE_POSITIVE(this->sum_of_mu_n_k);
		this->sum_of_mu_n_k += new_val;
	}
	void verify(const Q &q) const {
		long double verify_sum = 0.0L;
		For(node, q.Q_) {
			For(cell, *node) {
				verify_sum += *cell;
			}
		}
		assert(VERYCLOSE(this->sum_of_mu_n_k , verify_sum));
		this->sum_of_mu_n_k = verify_sum;
	}
};

struct Q_one_node : public Q :: Q_listener {
	mutable vector<long double> each_node;
	Q_one_node(const int N) : each_node(N) {}
	virtual void notify(int i, int, long double old_val, long double new_val, const Q *) {
		this->each_node.at(i) -= old_val;
		SHOULD_BE_POSITIVE(this->each_node.at(i));
		this->each_node.at(i) += new_val;
		if(this->each_node.at(i)>1.0L) {
			assert(VERYCLOSE(1, this->each_node.at(i)));
			this->each_node.at(i) = 1.0L;
		}
	}
	void verify(const Q &q) const {
		vector<long double> verify_each_node(q.N);
		for(int i=0; i<q.N; ++i) {
			For(cell, q.Q_.at(i)) {
				verify_each_node.at(i) += *cell;
			}
		}
		for(int i=0; i<q.N; ++i) {
			assert(VERYCLOSE(this->each_node.at(i) , verify_each_node.at(i)));
		}
		this->each_node = verify_each_node;
	}
};

namespace std {
namespace tr1 {

template<typename a, typename b>
struct hash< std::pair<a, b> > {
private:
   const hash<a> ah;
   const hash<b> bh;
public:
   hash() : ah(), bh() {}
   size_t operator()(const std::pair<a, b> &p) const {
      return ah(1+p.first) ^ bh(-p.second);
   }
};

}} // namespaces;
template struct std :: tr1 :: hash< pair<int,int> >;

template<int Power>
struct Q_templated_y_kl : public Q :: Q_listener {
	// Summing ONLY across the edges,
	// this will store Q_ik \times Q_jl,
	// or their square.
	// Must take care of direction too.
	// If directed, y_kl is the edges from k to l
	// If undirected, make k the smaller cluster-id
	VVL y_kl;
	Network * const net;
	Q_templated_y_kl(Network *net_) : y_kl(J, vector<long double>(J, 0.0L)), net(net_) {}

	virtual void notify(int i, int k, long double old_val, long double new_val, const Q *q) {
		// Must consider all the Junctions at node i
		// .. but beware of self loops.
		vector<long double> sum__of__Q_l__for__the__nonselfloop__out__neighbours__of__i(J);
		vector<long double> sum__of__Q_l__for__the__nonselfloop__in___neighbours__of__i(J);
		For(junc_id, net->i.at(i).my_junctions) {
			const Junction & junc = net->junctions->all_junctions_sorted.at(*junc_id);
			assert(junc.this_node_id == i);
			const int j = junc.far_node_id;
			if(j==i)
				continue; // we'll ignore self loops in this optimization for now
			for(int l = 0; l < J; ++l) {
				const long double Qjl = q->get(j,l);
				if(junc.junction_type == -1)
					sum__of__Q_l__for__the__nonselfloop__in___neighbours__of__i.at(l) += power<Power>(Qjl); // in
				else
					sum__of__Q_l__for__the__nonselfloop__out__neighbours__of__i.at(l) += power<Power>(Qjl); // out and undir
			}
		}
		for(int l = 0; l < J; ++l) {
			const long double total_out = sum__of__Q_l__for__the__nonselfloop__out__neighbours__of__i.at(l);
			const long double total_in  = sum__of__Q_l__for__the__nonselfloop__in___neighbours__of__i.at(l);
			if(net->directed) {
				this->y_kl.at(k).at(l) += (power<Power>(new_val) - power<Power>(old_val)) * total_out;
				this->y_kl.at(l).at(k) += (power<Power>(new_val) - power<Power>(old_val)) * total_in;
			} else { // undirected
				int k2 = k, l2 = l;
				if(k2>l2)
					swap(k2,l2);
				this->y_kl.at(k2).at(l2) += (power<Power>(new_val) - power<Power>(old_val)) * total_out;
			}
		}
		if(net->has_self_loop.at(i))
		{
			const int j = i;
			// each of these junctions corresponds to an edge
			// we should consider all possible clusters for the other endpoint
			for(int l = 0; l < J; ++l) {
				const long double Qjl = q->get(j,l);
				int k2 = k;
				int l2 = l;
				if(this->net->directed==false) { // if undirected, we should have k2 <= j2
					if(k2>l2)
						swap(k2,l2);
				}
				assert(j==i);
				{
					if(k==l) // y_kk
						this->y_kl.at(k2).at(l2) += power<Power>(new_val * new_val) - power<Power>(old_val * old_val);
					else { // y_kl, where k != l
						if(this->net->directed) {
							// we only deal with *one* of the two self-loop junctions,
							// therefore we have to double up here:
							this->y_kl.at(k2).at(l2) += power<Power>(new_val * Qjl) - power<Power>(old_val * Qjl);
							this->y_kl.at(l2).at(k2) += power<Power>(new_val * Qjl) - power<Power>(old_val * Qjl);
						} else {
							this->y_kl.at(k2).at(l2) += 2*power<Power>(new_val * Qjl) - 2*power<Power>(old_val * Qjl);
						}
					}
				}
			}
		}
	}
	void verify(const Q &q) const {
		VVL verify_y_kl(J, vector<long double>(J) );
		const int E = net->edge_set->E();
		for(int m=0; m<E; ++m) {
			EdgeSet :: Edge edge = net->edge_set->edges.at(m);
			int i = edge.left;
			int j = edge.right;
			for(int k=0; k<J; ++k) {
				for(int l=0; l<J; ++l) {
					const long double Qik = q.get(i,k);
					const long double Qjl = q.get(j,l);
					int k2 = k;
					int l2 = l;
					if(this->net->directed == false) { // if undirected, we should have k2 <= j2
						if(k2>l2)
							swap(k2,l2);
					}
					verify_y_kl.at(k2).at(l2) += power<Power>(Qik * Qjl);
				}
			}
		}
		// PP( this->y_kl.size() );
		assert(this->y_kl.size() == J);
		assert(verify_y_kl.size() == J);
		for(int k=0; k<J; ++k) {
			assert(this-> y_kl.at(k).size() == J);
			assert(verify_y_kl.at(k).size() == J);
			for(int l=0; l<J; ++l) {
				DYINGWORDS(VERYCLOSE(  this->y_kl.at(k).at(l)
			                  ,       verify_y_kl.at(k).at(l)
				)) {
					PP(this->y_kl.at(k).at(l) -       verify_y_kl.at(k).at(l) );
					PP5(Power, k,l,  this->y_kl.at(k).at(l) ,       verify_y_kl.at(k).at(l) );
				}
			}
		}
	}
};
typedef Q_templated_y_kl<1> Q_mu_y_kl;
typedef Q_templated_y_kl<2> Q_squared_y_kl;

template<int Power>
struct Q_templated_psl_kl : public Q :: Q_listener {
	VVL psl_kl;
	Q_templated_psl_kl(): psl_kl(J, vector<long double>(J,0)) {}

	virtual void notify(int i, int k, long double old_val, long double new_val, const Q *q) {
		for(int l=0; l<J; ++l) {
			if(l==k) {
				this->psl_kl.at(k).at(l) += power<Power>(new_val * new_val) - power<Power>(old_val * old_val);
			} else {
				const long double Qil = q->get(i,l);
				this->psl_kl.at(k).at(l) += power<Power>(new_val * Qil) - power<Power>(old_val * Qil);
				this->psl_kl.at(l).at(k) += power<Power>(new_val * Qil) - power<Power>(old_val * Qil);
			}
		}
	}
	void verify(const Q &q) const {
		VVL verify_psl_kl(J, vector<long double>(J,0));
		for(int i=0; i<q.N; ++i) {
			for(int k=0; k<J; ++k) {
				for(int l=0; l<J; ++l) {
					const long double Qik = q.get(i,k);
					const long double Qil = q.get(i,l);
					verify_psl_kl.at(k).at(l) += power<Power>(Qik * Qil);
				}
			}
		}
		assert(this->psl_kl.size() == J);
		assert(verify_psl_kl.size() == J);
		for(int k=0; k<J; ++k) {
			assert(this-> psl_kl.at(k).size() == J);
			assert(verify_psl_kl.at(k).size() == J);
			for(int l=0; l<J; ++l) {
				assert(VERYCLOSE(  this->psl_kl.at(k).at(l)
			                  ,       verify_psl_kl.at(k).at(l)
				));
			}
		}
	}
};
typedef Q_templated_psl_kl<1> Q_mu_psl_kl;
typedef Q_templated_psl_kl<2> Q_sq_psl_kl;

gsl_rng * global_r = NULL;

// A few globals (sorry) but they only help in verification and assertions
struct Tracker {
	const Q * const q;
	Network * const net;
	const Q_mu_n_k        * const ql_mu_n_k;
	const Q_squared_n_k   * const ql_squared_n_k;
	const Q_entropy       * const ql_entropy;
	const Q_sum_of_mu_n_k * const ql_sum_of_mu_n_k;
	const Q_one_node      * const ql_one_node     ;
	const Q_mu_y_kl       * const ql_mu_y_kl;
	const Q_squared_y_kl  * const ql_squared_y_kl;
	const Q_mu_psl_kl     * const ql_mu_psl_kl;
	const Q_sq_psl_kl     * const ql_sq_psl_kl;
	Tracker(const Q *q_, Network * net_
			, const Q_mu_n_k *ql_mu_n_k_
			, const Q_squared_n_k *ql_squared_n_k_
			, const Q_entropy *ql_entropy_
			, const Q_sum_of_mu_n_k *ql_sum_of_mu_n_k_
			, const Q_one_node *ql_one_node_
			, const Q_mu_y_kl *ql_mu_y_kl_
			, const Q_squared_y_kl *ql_squared_y_kl_
			, const Q_mu_psl_kl *ql_mu_psl_kl_
			, const Q_sq_psl_kl *ql_sq_psl_kl_
			)
		: q(q_), net(net_)
		  , ql_mu_n_k(ql_mu_n_k_)
		  , ql_squared_n_k(ql_squared_n_k_)
		  , ql_entropy(ql_entropy_)
		  , ql_sum_of_mu_n_k(ql_sum_of_mu_n_k_)
		  , ql_one_node(ql_one_node_)
		  , ql_mu_y_kl(ql_mu_y_kl_)
		  , ql_squared_y_kl(ql_squared_y_kl_)
		  , ql_mu_psl_kl(ql_mu_psl_kl_)
		  , ql_sq_psl_kl(ql_sq_psl_kl_)
	{}
	void verify_all() const {
		ql_mu_n_k->verify(*this->q);
		ql_squared_n_k->verify(*this->q);
		ql_entropy->verify(*this->q);
		ql_sum_of_mu_n_k->verify(*this->q);
		ql_one_node->verify(*this->q);
		ql_mu_y_kl->verify(*this->q);
		ql_squared_y_kl->verify(*this->q);
		ql_mu_psl_kl->verify(*this->q);
		ql_sq_psl_kl->verify(*this->q);
	}
};
const Tracker * global_tracker = NULL;

// Hyperparameters
const long double alpha_for_stick_breaking = 1.0L;
const long double beta_1 = 1.0L;
const long double beta_2 = 1.0L;


static long double exp_log_Taylor_Normal(const long double nonRandom, const long double mean_, /* MUTABLE */ long double variance);
static long double exp_log_Gamma_logNormal(const long double nonRandom, const long double mean, const long double variance);
#define exp_log_Gamma_Normal exp_log_Gamma_logNormal
// #define exp_log_Gamma_Normal exp_log_Taylor_Normal

struct SetOfNormalPercentiles {
	static const int num_percs = 10;
	array<long double, num_percs> perc;
	SetOfNormalPercentiles() {
		for(int i=0; i<num_percs; ++i) {
			const long double P = (1.0L+i*2.0L)/(num_percs*2.0L);
			perc.at(i) = gsl_cdf_gaussian_Pinv(P, 1);
			// PP3(i,P,perc.at(i));
		}
	}
};
static long double exp_log_Gamma_logNormal(const long double nonRandom, const long double mean, const long double variance) {
	if(variance == 0) {
		return gsl_sf_lngamma(mean + nonRandom);
	}
	assert(variance>0);
	static SetOfNormalPercentiles set_of_normal_percentiles;
	// First, identify the logNormal parameters, m and v
	const long double s2 = logl(variance / mean /mean + 1);
	const long double m = logl(mean) - s2/2.0L;
	const long double s = sqrtl(s2);
	assert(isfinite(s2));
	assert(isfinite(m));
	assert(isfinite(s));
	long double answer = 0.0L;
	// array<double, SetOfNormalPercentiles :: num_percs> log_normally_distributed;
	for(int i=0; i<set_of_normal_percentiles.num_percs; ++i) {
		const long double x = set_of_normal_percentiles.perc[i] * s + m;
		const long double expx = expl(x);
		// PP6(i,mean, sqrtl(variance), set_of_normal_percentiles.perc.at(i), x, expx);
		// log_normally_distributed.at(i) = expx;
		answer += gsl_sf_lngamma(nonRandom + expx);
	}
	// const double observed_mean = gsl_stats_mean(log_normally_distributed.data(), 1, SetOfNormalPercentiles :: num_percs);
	// const double observed_sd   = gsl_stats_sd  (log_normally_distributed.data(), 1, SetOfNormalPercentiles :: num_percs);
	// PP4(mean, observed_mean, sqrtl(variance), observed_sd);
	//For(lnd , log_normally_distributed) {
	//}
	answer /= set_of_normal_percentiles.num_percs;
	assert(isfinite(answer));
	// const long double old = exp_log_Taylor_Normal(nonRandom, mean, variance);
	// PP2(old, answer);
	// return old;
	return answer;
}
static long double exp_log_Taylor_Normal(const long double nonRandom, const long double mean_, /* MUTABLE */ long double variance) {
	// we're trying to calculcate E [ log ( Gamma ( nonRandom + N(mean,variance) ) ) ]
	if(variance > mean_) {
		assert(VERYCLOSE(variance, mean_));
		variance = mean_;
	}
	assert(variance >= 0.0L);
	assert(variance <= mean_);
	assert(mean_ >= 0.0L);
	if(mean_ == 0.0L) {
		assert(variance == 0.0L);
	}
	else {
		// we should ensure that the std dev is not too large compared to the mean.
		// In particular, I'm referring to the mean of the random bit
		const long double stddev = sqrtl(variance);
		const long double x = stddev / mean_;
		const long double new_x = x / (x+1.0);
		assert(new_x >= 0.0L && new_x <= 1.0L && new_x <= x);
		const long double newstddev = new_x * mean_;
		assert(newstddev <= stddev);
		const long double newvariance = newstddev * newstddev;
		assert(isfinite(newvariance));
		variance = newvariance;
	}
	const long double mean = nonRandom + mean_;
	static const long double half_logl2pi = 0.5L * logl(2 * M_PI);
	// assert(variance >= 0.0L);
	// assert(mean >= 0.0L);
	return (mean - 0.5L) * logl(mean)
		+ (0.5L) * variance * ( 1.0L/mean + 0.25L/mean/mean )
		- mean
		+ half_logl2pi;
}
static long double gamma_k(const int k) {
	// k = J-k -1;
	assert(k>=0);
	assert(k<J);
	if(k<FIXED_K) return 1; else return 0.001;
	return powl(alpha_for_stick_breaking / (1.0L+alpha_for_stick_breaking), k);
}

template<typename T>
long double mapat (const T &container, int k, int l) {
	const typename T :: const_iterator location = container.find(make_pair(k,l));
	if(location == container.end())
		return 0.0L;
	else
		return location->second;
}

static void dump_block_summary(bool known_full = false) {
	const int N = global_tracker->q->N;
	if(VERYCLOSE(N, global_tracker->ql_sum_of_mu_n_k->sum_of_mu_n_k))
		known_full = true;
	if(known_full)
		assert(VERYCLOSE(N, global_tracker->ql_sum_of_mu_n_k->sum_of_mu_n_k));
	const VVL &mu_y_kl = global_tracker->ql_mu_y_kl->y_kl;
	const vector<long double> & mu_n_k = global_tracker->ql_mu_n_k->n_k;
	enum X { Ratio=0, Percentage=1, Done=2 };
	cout << stack.push << fixed << setprecision(1);
for(int x = Percentage; x<Done; ++x) {
	cout << (x==Ratio?"            ":"            ");
	for(int l=0; l<J; ++l) {
		cout << setw(x==Ratio?12:6) << mu_n_k.at(l);
	}
	cout << endl;
	long double count_all_edges = 0.0L;
	long double count_all_pairs = 0.0L;
	for(int k=0; k<J; ++k) {
		cout << "<" << setw(3) << k << " > " << setw(5) << mu_n_k.at(k);
		for(int l=0; l<J; ++l) {
			int k2 = k;
			int l2 = l;
			if(!global_tracker->net->directed && k2>l2) {
				swap(k2,l2);
			}
			const long double mu_y_kl_ = mu_y_kl.at(k2).at(l2);
			long double mu_p_kl = mu_n_k.at(k2)*mu_n_k.at(l2);
			if(global_tracker->net->directed ==false) {
				mu_p_kl  += global_tracker->ql_mu_psl_kl->psl_kl.at(k).at(l);
			}
			if(global_tracker->net->directed == false && l2==k2) { // if undirected
				mu_p_kl /= 2;
			}
			long double nonEdge_mu = mu_p_kl - mu_y_kl_;
			SHOULD_BE_POSITIVE(nonEdge_mu);
			switch(x) {
				break; case Ratio: cout << '|' << setw(5) << mu_y_kl_ << '/' << setw(5) << mu_p_kl;
				break; case Percentage:
					if(VERYCLOSE(mu_p_kl,0.0L))
						cout << '|' << setw(5) << ' ';
					else
						cout << '|' << setw(5) << mu_y_kl_/mu_p_kl*100.0L;
			}

			if(global_tracker->net->directed == false) {
				if(k<=l) {
					count_all_edges += mu_y_kl_;
					count_all_pairs += mu_p_kl;
				}
			} else {
				count_all_edges += mu_y_kl_;
				count_all_pairs += mu_p_kl;
			}
		}
		cout << endl;
	}
	if(known_full) {
		const int N = global_tracker->q->N;
		const int E = global_tracker->net->E();
		const int expected_pairs = global_tracker->net->directed ? N*N : N*(N+1)/2;
		assert(VERYCLOSE(count_all_edges, E));
		if(global_tracker->net->directed)
			assert(VERYCLOSE(count_all_pairs, N*N));
		else
			assert(VERYCLOSE(count_all_pairs, N*(N+1)/2));
	}
}
	cout << stack.pop;
}

long double calculate_first_four_terms_slowly(const Q *
// #define SLOW_P_KL
#ifdef SLOW_P_KL
		q
#endif
		, Network *net
		, const bool verbose = false
		) {
	assert(global_tracker);
	// global_tracker->verify_all();
	const vector<long double> & mu_n_k = global_tracker->ql_mu_n_k->n_k;
	const vector<long double> & sq_n_k = global_tracker->ql_squared_n_k->n_k;

	const VVL &mu_y_kl = global_tracker->ql_mu_y_kl->y_kl;
	const VVL &sq_y_kl = global_tracker->ql_squared_y_kl->y_kl;

#ifdef SLOW_P_KL
	// These next structures should be removed sometimes.  For now,
	// they are just for paranoid testing of p_kl
	vector< vector<long double> > mu_slowp_kl(J, vector<long double>(J) );
	vector< vector<long double> > sq_slowp_kl(J, vector<long double>(J) );
	const int N = q->N;
	for(int i=0; i<N; ++i) {
	for(int j= (net->directed?0:i) ; j<N; ++j) { // undirected j starts at i
		for(int k=0; k<J; ++k) {
			for(int l=0; l<J; ++l) {
				const long double Qik = q->get(i,k);
				const long double Qjl = q->get(j,l);
				int k2 = k;
				int l2 = l;
				if(net->directed == false) { // if undirected, we should have k2 <= j2
					if(k2>l2)
						swap(k2,l2);
				}
				mu_slowp_kl.at(k2).at(l2) += Qik * Qjl;
				sq_slowp_kl.at(k2).at(l2) += Qik * Qjl * Qik * Qjl;
			}
		}
	}
	}

if(net->directed == false)
	for(int k=0; k<J; ++k) {
		for(int l=0; l<J; ++l) {
			if(l<k) {
				assert(0==mu_slowp_kl.at(k).at(l));
				assert(0==sq_slowp_kl.at(k).at(l));
				continue;
			}
			// PP2(mu_slowp_kl.at(k).at(l), sq_slowp_kl.at(k).at(l));
		}
	}
#endif

	long double four_terms_1cluster_sizes = 0.0L;
	long double four_terms_2edges = 0.0L;
	long double four_terms_3non_edges = 0.0L;
	long double four_terms_4pairs = 0.0L;

	// First term, E log Gamma( n_k + gamma_k )
	for(int k=0; k<J; k++) {
		long double mu = mu_n_k.at(k);
		long double var = mu - sq_n_k.at(k);
		SHOULD_BE_POSITIVE(mu);
		SHOULD_BE_POSITIVE(var);
		four_terms_1cluster_sizes += exp_log_Gamma_Normal( gamma_k(k), mu, var );
	}
	if(verbose)
		PP(four_terms_1cluster_sizes);
	// Second third and fourth terms, E log Gamma (y_kl + Beta_1)
//*
	// cout << endl;
	for(int k=0; k<J; k++) {
		for(int l=0; l<J; l++) {
			if(l<k && net->directed == false) {
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

			// if undirected, the self loops sort of count twice
			if(net->directed ==false) {
				mu_p_kl  += global_tracker->ql_mu_psl_kl->psl_kl.at(k).at(l);
				var_p_kl += global_tracker->ql_mu_psl_kl->psl_kl.at(k).at(l) - global_tracker->ql_sq_psl_kl->psl_kl.at(k).at(l);
			}

			if(net->directed == false && l==k) { // if undirected, the diagonal blocks should be halved
				mu_p_kl /= 2;
				var_p_kl /= 2;
			}
#ifdef SLOW_P_KL
			const long double mu_slowp_KL = mu_slowp_kl.at(k).at(l);
			const long double var_slowp_KL = mu_slowp_kl.at(k).at(l) - sq_slowp_kl.at(k).at(l);
#endif

			//for the non-edges
			long double nonEdge_mu = mu_p_kl - mu;
			long double nonEdge_var = var_p_kl - var_y_kl;

			/*
			cout << k << ',' << l
				<< '\t' << mu << '(' << var_y_kl << ')'
				<< '\t' << mu_p_kl << '(' << var_p_kl << ')'
#ifdef SLOW_P_KL
				<< '=' << mu_slowp_KL << '(' << var_slowp_KL << ')'
#endif
				<< '\t' << nonEdge_mu << '(' << nonEdge_var << ')'
				<< endl;
				// */

			SHOULD_BE_POSITIVE(mu);
			SHOULD_BE_POSITIVE(var_y_kl);
			SHOULD_BE_POSITIVE(mu_p_kl);
			SHOULD_BE_POSITIVE(var_p_kl);
			SHOULD_BE_POSITIVE(nonEdge_mu);
			SHOULD_BE_POSITIVE(nonEdge_var);
#ifdef SLOW_P_KL
			assert(VERYCLOSE(mu_slowp_KL , mu_p_kl));
			assert(VERYCLOSE(var_slowp_KL , var_p_kl));
#endif

			four_terms_2edges += exp_log_Gamma_Normal( beta_1, mu, var_y_kl );
			four_terms_3non_edges += exp_log_Gamma_Normal( beta_2, nonEdge_mu, nonEdge_var );
			four_terms_4pairs -= exp_log_Gamma_Normal( beta_1 + beta_2, mu_p_kl, var_p_kl );
		}
	}
// */

	// NOTE, we DO NOT include the entropy term in the return.

	const long double sum_four_terms = four_terms_1cluster_sizes + four_terms_2edges + four_terms_3non_edges + four_terms_4pairs;
	if(verbose) {
		PP3(sum_four_terms , four_terms_1cluster_sizes ,four_terms_2edges+four_terms_3non_edges+four_terms_4pairs);
		PP3(four_terms_2edges , four_terms_3non_edges , four_terms_4pairs);
	}
	return sum_four_terms;
}
long double lower_bound(const Q &q, Network *net) {
	return global_tracker->ql_entropy->entropy
		+ calculate_first_four_terms_slowly(&q, net);
}
// The code above calculates stuff, below we have the actual algorithm.

void one_node_all_k(Q *q, Network * net, const int node_id);

void one_random_node_all_k(Q *q, Network * net) {
	// pick one node at random, and call one_node_all_k on it
	const int N = q->N;
	const double unif = gsl_rng_uniform(global_r);
	const int random_node = N * unif;
	one_node_all_k(q, net, random_node);
}

static void vacate_a_node(Q *q, const int node_id);
static const vector<long double> vacate_a_node_and_calculate_its_scores(Q *q, Network *net, const int node_id, const vector<int> * some_clusters = NULL, const bool verbose = false);
static bool check_total_score_is_1(const vector<long double> &scores);

void one_node_all_k(Q *q, Network * net, const int node_id) {
	const vector<long double> scores = vacate_a_node_and_calculate_its_scores(q, net, node_id);
	assert((int)scores.size() == J);
	assert(check_total_score_is_1(scores));

	for(int k=0; k<J; ++k) {
		q->set(node_id, k) = scores.at(k);
	}
}
void one_node_all_k_M3(Q *q, Network * net, const int node_id, const vector<int> * some_clusters = NULL) {
	const bool verbose = false;
if(verbose) {
	if(global_groundTruth_ptr)
		PP(global_groundTruth.at(node_id));
	if(some_clusters) {
		For(cl, *some_clusters) {
			PP(*cl);
		}
	}
}
	const int num_clusters_to_consider = (some_clusters == NULL) ? J : some_clusters->size();
	// assign a node totally to one cluster selected at random
	const vector<long double> scores = vacate_a_node_and_calculate_its_scores(q, net, node_id, some_clusters, verbose);
	for(size_t l = 0; l < scores.size(); ++l) {
		const long double score = scores.at(l);
		const int cluster = some_clusters ? some_clusters->at(l) : l;
if(verbose)
		PP3(cluster, node_id, score);
	}
	assert((int)scores.size() == num_clusters_to_consider);
	assert(check_total_score_is_1(scores));

	long double unif = gsl_rng_uniform(global_r);
	int random_offset = 0;
	while(random_offset+1 < (int)scores.size() && unif > scores.at(random_offset)) {
		unif -= scores.at(random_offset);
		++ random_offset;
	}
	assert(random_offset>=0 && random_offset < num_clusters_to_consider);

	int new_cluster = -1;
	if(some_clusters == NULL) {
		new_cluster = random_offset;
	} else {
		new_cluster = some_clusters->at(random_offset);
	}
if(verbose)
	PP(new_cluster);
else {
	cout << ", " << new_cluster;
}
	q->set(node_id, new_cluster) = 1;
}

static bool check_total_score_is_1(const vector<long double> &scores) {
	long double check_new_total_is_1 = 0.0L;
	For(score, scores) {
		// PP(*score);
		check_new_total_is_1 += *score;
	}
	return VERYCLOSE(check_new_total_is_1, 1.0L);
}

int my_primary_cluster(const int node_id, const Q & q) {
	const vector<long double> & my_clusters = q.Q_.at(node_id);
	return max_element(my_clusters.begin(), my_clusters.end()) - my_clusters.begin();
}
vector<int> pick_random_clusters(const int num_clusters_with_replacement, const Q & q) {
	assert(num_clusters_with_replacement > 0);
	// pick a node at random, and use its cluster
	const int N = q.N;
	vector<int> nodes;
	for(int attempt = 0; attempt < num_clusters_with_replacement; ++attempt) {
		const int random_node = gsl_rng_uniform(global_r) * N;
		const int my_cluster = my_primary_cluster(random_node, q);
		assert(my_cluster >= 0 && my_cluster < J);
		for(int i=0; i<N; ++i) {
			const long double q_ik = q.get(i,my_cluster);
			if(
					q_ik > 0.3 ||
					gsl_rng_uniform(global_r) < q_ik)
				nodes.push_back(i);
		}
	}
	sort(nodes.begin(), nodes.end());
	vector<int> unique_nodes;
	unique_copy(nodes.begin(), nodes.end(), back_inserter(unique_nodes));
	random_shuffle(unique_nodes.begin(), unique_nodes.end());

	assert(unique_nodes.size() <= q.N);
	return unique_nodes;
}

static const vector<long double> vacate_a_node_and_calculate_its_scores(Q *q, Network *net, const int node_id, const vector<int> * some_clusters /*= NULL*/, const bool verbose /* = false*/) {
	vacate_a_node(q, node_id);
	// store the baseline ?
	const int num_clusters_to_consider = (some_clusters == NULL) ? J : some_clusters->size();
	vector<long double> scores(num_clusters_to_consider);
	for (int x = 0; x < num_clusters_to_consider; ++x) {
		const int target_cluster = (some_clusters == NULL) ? x : some_clusters->at(x);
		// cout << "trying node " << node_id << " in cluster " << k << endl;
		// assert(VERYCLOSE(0, global_tracker->ql_one_node->each_node.at(node_id)));
		q->set(node_id, target_cluster) = 1;
		// assert(VERYCLOSE(1, global_tracker->ql_one_node->each_node.at(node_id)));
		scores.at(x) = calculate_first_four_terms_slowly(q, net, verbose);
		//PP(scores.at(k));
		q->set(node_id, target_cluster) = 0;
		// exit(1);
	}
	assert(0 == global_tracker->ql_one_node->each_node.at(node_id));
	const long double max_score = *max_element(scores.begin(), scores.end());
	for(int k=0; k< (int)scores.size(); ++k) {
		long double & score = scores.at(k);
		// PP2(k, score);
		score -= max_score;
		// const long double jitter=gsl_ran_gaussian(global_r, 0.0001); *score += jitter;
	}
	long double total = 0.0L;
	For(score, scores) {
		total += expl(*score);
	}
	For(score, scores) {
		*score  = expl(*score) / total;
	}
	// assert(0 == global_tracker->ql_one_node->each_node.at(node_id));
	return scores;
}
static void vacate_a_node(Q *q, const int node_id) {
	const int num_clusters = q->Q_.at(node_id).size();
	assert(J==num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		q->set(node_id, k) = 0;
	}
}

static vector<int> random_list_of_all_nodes(const int N) {
	vector<int> all_nodes_randomly;
	for(int i=0; i<N; i++)
		all_nodes_randomly.push_back(i);
	random_shuffle(all_nodes_randomly.begin(), all_nodes_randomly.end());
	return all_nodes_randomly;
}

static void Var_on_all_nodes(Q *q, Network *net, const int repetitions = 1) {
	vector<int> all_nodes_randomly = random_list_of_all_nodes(q->N);
	for(int rep = 0; rep < repetitions; ++rep) {
		For(i, all_nodes_randomly) {
			one_node_all_k(q, net, *i);
		}
	}
}

static void vacate_everything_then_M3_then_a_few_Var_moves(Q *q, Network * net) {
	const int N = q->N;
	for(int i=0; i<N; i++) {
		vacate_a_node(q, i);
	}
	assert(VERYCLOSE(0, global_tracker->ql_sum_of_mu_n_k->sum_of_mu_n_k));
	vector<int> all_nodes_randomly = random_list_of_all_nodes(N);
	For(i, all_nodes_randomly) {
		one_node_all_k_M3(q, net, *i);
	}
	for(int multiVar=0; multiVar < 20; ++multiVar) {
		vector<int> all_nodes_randomly = random_list_of_all_nodes(N);
		For(i, all_nodes_randomly) {
			one_node_all_k(q, net, *i);
		}
	}
	// dump(q,net);
	dump_block_summary();
	global_tracker->ql_mu_n_k->dump_me();
}
static void vacate_somenodes_then_M3_then_a_few_Var_moves(Q *q, Network * net, const vector<int> &random_nodes, int howManyVar) {
	cout << "vacating... with " << random_nodes.size() << " nodes." << endl;
	const int N = q->N;
	For(i, random_nodes) {
		vacate_a_node(q, *i);
	}
	// cout << "should be some missing now" << endl; dump_block_summary();
	dump_block_summary(false); cout << "  vacated" << endl;
	For(i, random_nodes) {
		one_node_all_k_M3(q, net, *i);
	}
	dump_block_summary(true); cout << "  M3" << endl;
	for(int multiVar=0; multiVar < howManyVar; ++multiVar) {
		For(i, random_nodes) {
			one_node_all_k(q, net, *i);
		}
	}
	// dump_block_summary(true); cout << "  Var" << endl;
	// dump(q,net);
	// dump_block_summary();
	// global_tracker->ql_mu_n_k->dump_me();
}

void empty_one_cluster_then_M3_all_nodes_then_Var_all_nodes(Q &q, Network *net) {
	const int N = q.N;
		// Pick one clusters randomly, and set q_*k = 0
		dump_block_summary(true);
		const int random_node = gsl_rng_uniform(global_r) * N;
		const int k = my_primary_cluster(random_node, q);
		long double amount_removed = 0.0L;
		for(int i=0; i<N; ++i) {
			amount_removed += q.get(i,k);
			q.set(i,k) = 0;
		}
		dump_block_summary(false);
		cout << "emptied one cluster (ALL nodes) " << amount_removed << endl;
		vector<int> all_nodes_randomly = random_list_of_all_nodes(N);
		For(i, all_nodes_randomly) {
			one_node_all_k_M3(&q, net, *i);
		}
		dump_block_summary(true);
		cout << "M3 (ALL nodes)" << endl;

		Var_on_all_nodes(&q, net, 5);
		dump_block_summary(true);
		cout << "all nodes Var(x5)" << endl;
}

int make_random_choice_from_scores(const vector<long double> &scores) {
	long double unif = gsl_rng_uniform(global_r);
	int random_cluster = 0;
	while(random_cluster+1 < (int)scores.size() && unif > scores.at(random_cluster)) {
		unif -= scores.at(random_cluster);
		++ random_cluster;
	}
	return random_cluster;
}

void discretize_all_nodes(Q &q, Network * net) {
	const int N = q.N;
	for(int i=0; i<N; ++i) {
		const vector<long double> &scores = q.Q_.at(i);
		const int random_cluster = make_random_choice_from_scores(scores);
		for(int k=0; k<J; ++k) {
			q.set(i,k) = 0;
		}
		q.set(i,random_cluster) = 1;
	}
}

void discretize_then_M3(Q &q, Network * net) {
	discretize_all_nodes(q, net);
	dump_block_summary(true);
	cout << "discretized all nodes" << endl;

	const int N = q.N;
	const int rand_node = N * gsl_rng_uniform(global_r);
	int rand_cluster = my_primary_cluster(rand_node, q);

	// first, find if there are any other empty clusters.  If not, we must abandon.
	int another_empty_cluster = 0;
	while(true) {
		const long double n_k = global_tracker->ql_mu_n_k->n_k.at(another_empty_cluster);
		if(n_k < 0.5) {
			assert(VERYCLOSE(n_k, 0.0L));
			if(another_empty_cluster != rand_cluster)
				break;
		}
		++another_empty_cluster;
	}

	if(another_empty_cluster < rand_cluster)
		swap(another_empty_cluster, rand_cluster);

	if(another_empty_cluster >= J) {
		cout << "abandoning after discretization, due to lack of empty clusters "; PP(lower_bound(q,net));
		return;
	}
	assert(another_empty_cluster < J);
	assert(another_empty_cluster > rand_cluster);

	// next, empty that cluster
	cout << "             " << setw(rand_cluster*6) << "" << "<<0>>" << endl;
	const int n_rand = global_tracker->ql_mu_n_k->n_k.at(rand_cluster) + 0.1; // add 0.1 so that the it rounds correctly
	vector<int> nodes_in_the_random_cluster;
	{
		int verify_n_rand = 0;
		for(int i=0; i<N; ++i) {
			const long double Qir = q.get(i,rand_cluster);
			assert(VERYCLOSE(Qir,0) || VERYCLOSE(Qir,1));
			if(VERYCLOSE(Qir,1)) {
				++ verify_n_rand;
				nodes_in_the_random_cluster.push_back(i);
			}
			q.set(i,rand_cluster) = 0;
		}
		dump_block_summary(false);
		cout << "emptied one cluster" << endl;
		assert(n_rand == verify_n_rand);
	}

	assert(VERYCLOSE(0.0L,global_tracker->ql_mu_n_k->n_k.at(another_empty_cluster)));
	assert(VERYCLOSE(0.0L,global_tracker->ql_mu_n_k->n_k.at(rand_cluster)));
	PP2(rand_cluster, another_empty_cluster);
	vector<int> two_empty_clusters;
	two_empty_clusters.push_back(rand_cluster);
	two_empty_clusters.push_back(another_empty_cluster);
	assert(2==two_empty_clusters.size());

	vector<int> all_nodes_randomly = random_list_of_all_nodes(N);
	random_shuffle(nodes_in_the_random_cluster.begin(), nodes_in_the_random_cluster.end());
	PP(nodes_in_the_random_cluster.size());
	// For(i, all_nodes_randomly)
	cout << "m3 on those two clusters: ";
	For(i, nodes_in_the_random_cluster)
	{
		one_node_all_k_M3(&q, net, *i, &two_empty_clusters);
	}
	cout << endl;
	dump_block_summary(true);
	cout << "all nodes M3ed. "; PP(lower_bound(q,net));
}
bool swap_big_clusters_to_the_left(Q&q, Network *net) {
			for(int k=0; k+1 < J; ++k) {
				if (global_tracker->ql_mu_n_k->n_k.at(k)+0.01 < global_tracker->ql_mu_n_k->n_k.at(k+1)) {
					// PP2(global_tracker->ql_mu_n_k->n_k.at(k), global_tracker->ql_mu_n_k->n_k.at(k+1));
					const long double backup_score = global_tracker->ql_entropy->entropy + calculate_first_four_terms_slowly(&q, net);
					cout << "Swap " << k << " and " << k+1 << "...  ";
					for(int i=0; i<q.N; ++i) {
						const long double Q_k_i = q.get(i,k);
						const long double Q_k2i = q.get(i,k+1);
						q.set(i,k) = 0;
						q.set(i,k+1) = 0;
						q.set(i,k) = Q_k2i;
						q.set(i,k+1) = Q_k_i;
					}
					const long double swapped_score = global_tracker->ql_entropy->entropy + calculate_first_four_terms_slowly(&q, net);
					// PP3(backup_score, swapped_score, backup_score-swapped_score);
					assert( backup_score <= swapped_score + 0.001 );
					return true; // lambda return;
				}
			}
			return false;
}

void vcsbm(Network * net) {
	const int N = net->N();

	global_r = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set(global_r, args_info.seed_arg);

	// test_exp_log_Gamma_Normal(); return;

	Q q(N,J);
	Q_mu_n_k ql_mu_n_k;
	Q_squared_n_k ql_squared_n_k;
	Q_entropy ql_entropy;
	Q_sum_of_mu_n_k ql_sum_of_mu_n_k;
	Q_one_node ql_one_node(N);
	Q_mu_y_kl ql_mu_y_kl(net);
	Q_squared_y_kl ql_squared_y_kl(net);
	Q_mu_psl_kl ql_mu_psl_kl;
	Q_sq_psl_kl ql_sq_psl_kl;

	q.add_listener(&ql_mu_n_k);
	q.add_listener(&ql_squared_n_k);
	q.add_listener(&ql_entropy);
	q.add_listener(&ql_sum_of_mu_n_k);
	q.add_listener(&ql_one_node);
	q.add_listener(&ql_mu_y_kl);
	q.add_listener(&ql_squared_y_kl);
	q.add_listener(&ql_mu_psl_kl);
	q.add_listener(&ql_sq_psl_kl);

	assert(global_tracker == NULL);
	global_tracker = new Tracker(&q, net, &ql_mu_n_k, &ql_squared_n_k, &ql_entropy, &ql_sum_of_mu_n_k, &ql_one_node
			, &ql_mu_y_kl, &ql_squared_y_kl
			, &ql_mu_psl_kl, &ql_sq_psl_kl
			);
	assert(global_tracker);

	// To store the best one found so far.
	Q q_copy(N,J);
	long double best_lower_bound_found = - numeric_limits<long double> :: max();

	ql_entropy.verify(q);
	ql_mu_n_k.dump_me();
	ql_squared_n_k.dump_me();
	ql_entropy.verify(q);
	ql_mu_n_k.verify(q);
	ql_squared_n_k.verify(q);

	if(1)
	for(int i=0; i<N; i++) {
		// if(i%2==0) q.set(i,0) = 1; else q.set(i,1) = 1;
		// if(i<6) q.set(i,0) = 1; else q.set(i,1) = 1;
		// if(i>0) q.set(i,(i-1)/4)=1;
		// if(i>0) q.set(i,(i-1)/4)=1;
		// q.set(i,gsl_rng_uniform(global_r)*3)=1;
		// q.set(i,0) = 0.5; q.set(i,1) = 0.5;
		// one_node_all_k(&q, net, i);
		q.set(i,0) = 1; // everything in the first cluster.
		// q.set(i,i/10) = 1; // everything in the correct cluster at first // alpha=1, K10O10 => -2109.42
		// q.set(i,i/10) = 1; // everything in the correct cluster at first // alpha=3, K10O10 => -2118.26
		// q.set(i,i/10) = 1; // everything in the correct cluster at first // alpha=10, K10O10 => -2119.85
	}
	// dump(&q, net);
	dump_block_summary();
	const long double initial_score = ql_entropy.entropy + calculate_first_four_terms_slowly(&q, net);
	cout << "That was the initial state "; PP(initial_score);
	cout << endl;
	for(int repeat=0; repeat < 4; ++repeat) {
		PP(repeat);
		// do some changes, but revert if things don't improve
		VVL q_backupQ_ = q.Q_;

		const long double backup_score = ql_entropy.entropy + calculate_first_four_terms_slowly(&q, net);
#if 0
		vector<int> some_random_clusters = pick_random_clusters(gsl_rng_uniform(global_r) < 0.5 ? 1 : 2, q);

		// merge/split one or two clusters
		vacate_somenodes_then_M3_then_a_few_Var_moves(&q, net, some_random_clusters, 3);
		// a sweep on all the nodes
		dump_block_summary(true);
		cout << "one or two clusters M3+Var" << endl;

		Var_on_all_nodes(&q, net);
		dump_block_summary(true);
		cout << "all nodes Var" << endl;
#endif
		while(swap_big_clusters_to_the_left(q, net)) {
		}
		cout << endl;


		// empty_one_cluster_then_M3_all_nodes_then_Var_all_nodes(q,net);
		if(drand48() < 0.0) {
			discretize_then_M3(q,net);
		} else {
			Var_on_all_nodes(&q, net, 3);
			cout << "a Var(x3) on all nodes. "; PP(lower_bound(q,net));
		}
		dump_block_summary(true);


//#if 0
		const long double new_score = ql_entropy.entropy + calculate_first_four_terms_slowly(&q, net);
		if(new_score < backup_score) {
			cout << "Undo. From " << new_score << " back to " << backup_score << endl;
			for(int i=0; i<N; ++i) {
				vacate_a_node(&q, i);
				for(int k=0; k<J; ++k) {
					q.set(i,k) = q_backupQ_.at(i).at(k);
				}
			}
			assert(VERYCLOSE(backup_score , ql_entropy.entropy + calculate_first_four_terms_slowly(&q, net)));
		} else {
			cout << "Improved lower bound (repeat=" << repeat << "). Lower bound is: " << new_score << " was: " << backup_score << endl;
		}
//#endif

		cout << endl;
	}
	global_tracker->verify_all();
}
