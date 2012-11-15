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
#include<limits>


#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace std :: tr1;
using namespace network;

void vcsbm(Network * net);

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

#define FIXED_K 4
const int J = FIXED_K; // fix the upper bound on K at 10.

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
		cout << "node_name\t";
		cout << "degree\t";
		for(int k=0; k<J; ++k) {
			cout << "\tn_" << k;
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
	unordered_map< pair<int,int> , long double> y_kl;
	Network * const net;
	Q_templated_y_kl(Network *net_) : net(net_) {}

	virtual void notify(int i, int k, long double old_val, long double new_val, const Q *q) {
		// Must consider all the Junctions at node i
		// .. but beware of self loops.
		For(junc_id, net->i.at(i).my_junctions) {
			const Junction & junc = net->junctions->all_junctions_sorted.at(*junc_id);
			assert(junc.this_node_id == i);
			// each of these junctions corresponds to an edge
			// we should consider all possible clusters for the other endpoint
			for(int l = 0; l < J; ++l) {
				const int j = junc.far_node_id;
				const long double Qjl = q->get(j,l);
				int k2 = k;
				int l2 = l;
				{ // if undirected, we should have k2 <= j2
					if(k2>l2)
						swap(k2,l2);
				}
				this->y_kl[make_pair(k2,l2)] += power<Power>(new_val * Qjl) - power<Power>(old_val * Qjl);
			}
		}
	}
	void verify(const Q &q) const {
		vector< vector<long double> > verify_y_kl(J, vector<long double>(J) );
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
					{ // if undirected, we should have k2 <= j2
						if(k2>l2)
							swap(k2,l2);
					}
					verify_y_kl.at(k2).at(l2) += power<Power>(Qik * Qjl);
				}
			}
		}
		// PP( this->y_kl.size() );
		For( entry, this->y_kl ) {
			// PP2(               entry->second   ,  verify_y_kl.at(entry->first.first).at(entry->first.second)  );
			assert(VERYCLOSE(  entry->second   ,  verify_y_kl.at(entry->first.first).at(entry->first.second) ));
			verify_y_kl.at(entry->first.first).at(entry->first.second) = 0.0L;
		}
		For(verify_row, verify_y_kl) {
			const long double sum = accumulate(verify_row->begin(), verify_row->end(), 0);
			assert(VERYCLOSE(sum, 0.0L));
		}
	}
};
typedef Q_templated_y_kl<1> Q_mu_y_kl;
typedef Q_templated_y_kl<2> Q_squared_y_kl;

template<int Power>
struct Q_templated_psl_kl : public Q :: Q_listener {
	unordered_map< pair<int,int> , long double> psl_kl;

	virtual void notify(int i, int k, long double old_val, long double new_val, const Q *q) {
		for(int l=0; l<J; ++l) {
			if(l==k) {
				this->psl_kl[make_pair(k,l)] += power<Power>(new_val * new_val) - power<Power>(old_val * old_val);
			} else {
				const long double Qil = q->get(i,l);
				this->psl_kl[make_pair(k,l)] += power<Power>(new_val * Qil) - power<Power>(old_val * Qil);
				this->psl_kl[make_pair(l,k)] += power<Power>(new_val * Qil) - power<Power>(old_val * Qil);
			}
		}
	}
	void verify(const Q &q) const {
		unordered_map< pair<int,int> , long double> verify_psl_kl;
		for(int i=0; i<q.N; ++i) {
			for(int k=0; k<J; ++k) {
				for(int l=0; l<J; ++l) {
					const long double Qik = q.get(i,k);
					const long double Qil = q.get(i,l);
					verify_psl_kl[make_pair(k,l)] += power<Power>(Qik * Qil);
				}
			}
		}
		For(cell, this->psl_kl) {
			const int k = cell->first.first;
			const int l = cell->first.second;
			const long double val = cell->second;
			assert(VERYCLOSE(verify_psl_kl[make_pair(k,l)] , val));
		}
		For(cell, verify_psl_kl) {
			const int k = cell->first.first;
			const int l = cell->first.second;
			const long double val = cell->second;
			assert(VERYCLOSE(mapat(this->psl_kl,k,l) , val));
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
	const Q_mu_y_kl       * const ql_mu_y_kl;
	const Q_squared_y_kl  * const ql_squared_y_kl;
	const Q_mu_psl_kl     * const ql_mu_psl_kl;
	const Q_sq_psl_kl     * const ql_sq_psl_kl;
	Tracker(const Q *q_, Network * net_
			, const Q_mu_n_k *ql_mu_n_k_
			, const Q_squared_n_k *ql_squared_n_k_
			, const Q_entropy *ql_entropy_
			, const Q_sum_of_mu_n_k *ql_sum_of_mu_n_k_
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
		ql_mu_y_kl->verify(*this->q);
		ql_squared_y_kl->verify(*this->q);
		ql_mu_psl_kl->verify(*this->q);
		ql_sq_psl_kl->verify(*this->q);
	}
};
Tracker * global_tracker = NULL;

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
static long double gamma_k(const int k) {
	assert(k>=0);
	assert(k<J);
	if(k<FIXED_K) return 1; else return 0.001;
	return powl(alpha_for_stick_breaking / (1.0L+alpha_for_stick_breaking), k);
}

struct BreakdownOfCompleteRecalculation {
	BreakdownOfCompleteRecalculation() {}
	BreakdownOfCompleteRecalculation & reset() {
		assert(global_tracker);
		return *this;
	}
	void test_assuming_full() const {
		assert(global_tracker);
		const int N = global_tracker->q->N;
		assert(VERYCLOSE(N, global_tracker->ql_sum_of_mu_n_k->sum_of_mu_n_k));
	}
};

template<typename T>
long double mapat (const T &container, int k, int l) {
	const typename T :: const_iterator location = container.find(make_pair(k,l));
	if(location == container.end())
		return 0.0L;
	else
		return location->second;
}

long double calculate_first_four_terms_slowly(const Q *
#ifdef SLOW_P_KL
		q
#endif
		, Network * , BreakdownOfCompleteRecalculation &breakdown) {
	assert(global_tracker);
	global_tracker->verify_all();
	breakdown.reset();
	const vector<long double> & mu_n_k = global_tracker->ql_mu_n_k->n_k;
	const vector<long double> & sq_n_k = global_tracker->ql_squared_n_k->n_k;

	const unordered_map< pair<int,int> , long double> &mu_y_kl = global_tracker->ql_mu_y_kl->y_kl;
	const unordered_map< pair<int,int> , long double> &sq_y_kl = global_tracker->ql_squared_y_kl->y_kl;

#ifdef SLOW_P_KL
	// These next structures should be removed sometimes.  For now,
	// they are just for paranoid testing of p_kl
	vector< vector<long double> > mu_slowp_kl(J, vector<long double>(J) );
	vector< vector<long double> > sq_slowp_kl(J, vector<long double>(J) );
	const int N = q->N;
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

	long double first_4_terms = 0.0L;

	// First term, E log Gamma( n_k + gamma_k )
	for(int k=0; k<J; k++) {
		const long double mu = mu_n_k.at(k);
		const long double var = mu - sq_n_k.at(k);
		assert(mu >= 0.0L);
		assert(var >= 0.0L);
		first_4_terms += exp_log_Gamma_Normal( mu + gamma_k(k), var );
	}
	// Second third and fourth terms, E log Gamma (y_kl + Beta_1)
//*
	// cout << endl;
	for(int k=0; k<J; k++) {
		for(int l=0; l<J; l++) {
			if(l<k) {
				// if undirected, these block don't really count.
				// In particular y_kl should be zero
				assert(0 == mapat(mu_y_kl,k,l));
				assert(0 == mapat(sq_y_kl,k,l));
				continue;
			}
			long double mu = mapat(mu_y_kl,k,l);
			long double var_y_kl = mu - mapat(sq_y_kl,k,l);

			long double mu_p_kl = mu_n_k.at(k)*mu_n_k.at(l);
			long double var_p_kl = mu_n_k.at(k)*mu_n_k.at(l) - sq_n_k.at(k) * sq_n_k.at(l);
			mu_p_kl  += mapat(global_tracker->ql_mu_psl_kl->psl_kl,k,l);
			var_p_kl += mapat(global_tracker->ql_mu_psl_kl->psl_kl,k,l) - mapat(global_tracker->ql_sq_psl_kl->psl_kl,k,l);
			if(l==k) { // if undirected
				// Add in the self loops and then divide everything by two
				mu_p_kl /= 2;
				var_p_kl /= 2;
			} else { // bizarrelly, the self loops can contribute
			         // to between-cluster blocks.
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

void one_node_all_k(Q *q, Network * net, const int node_id);

void one_random_node_all_k(Q *q, Network * net) {
	// pick one node at random, and call one_node_all_k on it
	const int N = q->N;
	const double unif = gsl_rng_uniform(global_r);
	const int random_node = N * unif;
	one_node_all_k(q, net, random_node);
}

static void vacate_a_node(Q *q, const int node_id);
static const vector<long double> vacate_a_node_and_calculate_its_scores(Q *q, Network *net, const int node_id);
static bool check_total_score_is_1(const vector<long double> &scores);

void one_node_all_k(Q *q, Network * net, const int node_id) {
	const vector<long double> scores = vacate_a_node_and_calculate_its_scores(q, net, node_id);
	assert((int)scores.size() == J);
	assert(check_total_score_is_1(scores));

	for(int k=0; k<J; ++k) {
		q->set(node_id, k) = scores.at(k);
	}
}
void one_node_all_k_M3(Q *q, Network * net, const int node_id) {
	// assign a node totally to one cluster selected at random
	const vector<long double> scores = vacate_a_node_and_calculate_its_scores(q, net, node_id);
	assert((int)scores.size() == J);
	assert(check_total_score_is_1(scores));

	long double unif = gsl_rng_uniform(global_r);
	int random_cluster = 0;
	while(random_cluster+1 < (int)scores.size() && unif > scores.at(random_cluster)) {
		unif -= scores.at(random_cluster);
		++ random_cluster;
	}

	PP2(node_id, random_cluster);
	q->set(node_id, random_cluster) = 1;
}

static bool check_total_score_is_1(const vector<long double> &scores) {
	long double check_new_total_is_1 = 0.0L;
	For(score, scores) {
		// PP(*score);
		check_new_total_is_1 += *score;
	}
	return VERYCLOSE(check_new_total_is_1, 1.0L);
}

static const vector<long double> vacate_a_node_and_calculate_its_scores(Q *q, Network *net, const int node_id) {
	vacate_a_node(q, node_id);
	BreakdownOfCompleteRecalculation breakdown;
	// store the baseline ?
	vector<long double> scores(J);
	const int num_clusters = q->Q_.at(node_id).size();
	assert(J==num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		// cout << "trying node " << node_id << " in cluster " << k << endl;
		q->set(node_id, k) = 1;
		scores.at(k) = calculate_first_four_terms_slowly(q, net, breakdown);
		//PP(scores.at(k));
		q->set(node_id, k) = 0;
		// exit(1);
	}
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
	return scores;
}
static void vacate_a_node(Q *q, const int node_id) {
	const int num_clusters = q->Q_.at(node_id).size();
	assert(J==num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		q->set(node_id, k) = 0;
	}
}

static void vacate_everything_then_M3_then_a_few_Var_moves(Q *q, Network * net) {
	global_tracker->verify_all();
	BreakdownOfCompleteRecalculation breakdown;
	const int N = q->N;
	for(int i=0; i<N; i++) {
		vacate_a_node(q, i);
	}
	assert(VERYCLOSE(0, global_tracker->ql_sum_of_mu_n_k->sum_of_mu_n_k));
	calculate_first_four_terms_slowly(q, net, breakdown);
	for(int i=0; i<N; i++) {
		one_node_all_k_M3(q, net, i);
	}
	dump(q,net);
	calculate_first_four_terms_slowly(q, net, breakdown); breakdown.test_assuming_full();
	for(int i=0; i<N; i++) {
		one_node_all_k(q, net, i);
	}
	dump(q,net);
	global_tracker->ql_mu_n_k->dump_me();
	calculate_first_four_terms_slowly(q, net, breakdown); breakdown.test_assuming_full();
	global_tracker->verify_all();
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
	Q_mu_y_kl ql_mu_y_kl(net);
	Q_squared_y_kl ql_squared_y_kl(net);
	Q_mu_psl_kl ql_mu_psl_kl;
	Q_sq_psl_kl ql_sq_psl_kl;

	q.add_listener(&ql_mu_n_k);
	q.add_listener(&ql_squared_n_k);
	q.add_listener(&ql_entropy);
	q.add_listener(&ql_sum_of_mu_n_k);
	q.add_listener(&ql_mu_y_kl);
	q.add_listener(&ql_squared_y_kl);
	q.add_listener(&ql_mu_psl_kl);
	q.add_listener(&ql_sq_psl_kl);

	assert(global_tracker == NULL);
	global_tracker = new Tracker(&q, net, &ql_mu_n_k, &ql_squared_n_k, &ql_entropy, &ql_sum_of_mu_n_k
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

	BreakdownOfCompleteRecalculation breakdown;

	if(0)
	for(int i=0; i<N; i++) {
		// if(i%2==0) q.set(i,0) = 1; else q.set(i,1) = 1;
		// if(i<6) q.set(i,0) = 1; else q.set(i,1) = 1;
		// if(i>0) q.set(i,(i-1)/4)=1;
		// if(i>0) q.set(i,(i-1)/4)=1;
		// q.set(i,gsl_rng_uniform(global_r)*3)=1;
		// q.set(i,0) = 0.5; q.set(i,1) = 0.5;
		// one_node_all_k(&q, net, i);
		ql_mu_n_k.dump_me();
	}
	dump(&q, net);
	// everything assigned somewhere
	// calculate_first_four_terms_slowly(&q, net, true);
for(int restart = 0; restart<3; ++restart) {
	PP(restart);
	global_tracker->verify_all();
	if(0) {
		cout << endl << "random (re)initialization" << endl;
		for(int i=0; i<N; i++) {
			double rand_unif1 = gsl_rng_uniform(global_r);
			double rand_unif2 = gsl_rng_uniform(global_r);
			if(rand_unif1 > rand_unif2)
				swap(rand_unif1, rand_unif2);
			q.set(i,0) = 0.0; q.set(i,1) = 0.0; q.set(i,2) = 0.0;
			q.set(i,0) = rand_unif1; q.set(i,1) = rand_unif2-rand_unif1; q.set(i,2) = 1-rand_unif2;
		}
		dump(&q, net);
	}
	cout << endl << endl << "into the repeats now" << endl;
	for(int repeat = 0; repeat < 100; ++repeat) {
		PP2(restart,repeat);
		vacate_everything_then_M3_then_a_few_Var_moves(&q, net);
		const long double lower_bound = ql_entropy.entropy + calculate_first_four_terms_slowly(&q, net, breakdown);
		PP(lower_bound);
		{
			if(best_lower_bound_found < lower_bound) {
				cout << "New best lower bound found" << endl;
				best_lower_bound_found = lower_bound;
				q_copy.Q_ = q.Q_;
				dump(&q, net);
				ql_mu_n_k.dump_me();
				PP3(best_lower_bound_found, restart, repeat);
			}
		}
	}
}
}
