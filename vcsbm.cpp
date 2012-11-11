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
struct Q_n_k : public Q :: Q_listener {
	vector<long double> n_k;
	Q_n_k() : n_k(10000) {}
	virtual void notify(int, int k, long double old_val, long double new_val) {
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
	virtual ~Q_n_k() {}
};

void vcsbm(const Network * net) {
	const int N = net->N();
	const int J = 10; // fix the upper bound on K at 10.
	Q q(N,J);
	Q_n_k n_k;
	q.add_listener(&n_k);
	q.set(0,3) = 0.3;
	PP(q.get(0,3));
	q.set(0,3) = 0.6;
	PP(q.get(0,3));
	n_k.dump_me();
}
