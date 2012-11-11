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

#define assert_0_to_1(x) do { assert((x)>=0.0L); assert((x)<=1.0L); } while(0)

using namespace std;
using namespace network;

void vcsbm(const Network * net);

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
		assert_0_to_1(this->n_k.at(k));
		this->n_k.at(k) += new_val;
		assert_0_to_1(this->n_k.at(k));
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
	long double entropy;
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
		assert(this->entropy <= 0.0L);
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
		PP2(verify_entropy , this->entropy);
		assert(verify_entropy == this->entropy);
	}
};

void vcsbm(const Network * net) {
	const int N = net->N();
	const int J = 10; // fix the upper bound on K at 10.

	Q q(N,J);
	Q_mu_n_k mu_n_k;
	Q_squared_n_k squared_n_k;
	Q_entropy entropy;

	PP(entropy.entropy);

	q.add_listener(&mu_n_k);
	q.add_listener(&squared_n_k);
	q.add_listener(&entropy);
	entropy.verify(q);
	q.set(0,3) = 0.3;
	PP(q.get(0,3));
	q.set(0,3) = 0.6;
	PP(q.get(0,3));
	mu_n_k.dump_me();
	squared_n_k.dump_me();
	PP(entropy.entropy);
	entropy.verify(q);
	mu_n_k.verify(q);
	squared_n_k.verify(q);
}
