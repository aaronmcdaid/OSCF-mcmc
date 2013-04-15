#include"onmi.hpp"

/*
 * This code won't make much sense unless you use the arXiv document alongside it.
 *              http://arxiv.org/abs/1110.2515
 * As much as possible, I've tried to use the same variable names as there.
 */

#include<tr1/unordered_map>
#include<cmath>
#include<limits>
using namespace std :: tr1;

#include"macros.hpp"
#include"lvalue_input.hpp"
#include"hash_a_pair.hpp"
using namespace lvalue_input;

template <typename T>
void ignore(const T&) {
}

using namespace std;

long double h(const size_t w, const size_t n) {
	if(w==0)
		return 0.0L;
	return - (w * (log2l(w) - log2l(n)));
}

long double H_star(const size_t a, const size_t b, const size_t c, const size_t d, const size_t n) {
	const long double han = h(a,n);
	const long double hbn = h(b,n);
	const long double hcn = h(c,n);
	const long double hdn = h(d,n);
	if(han + hdn >= hbn + hcn) { // This is a normal, positively-correlated, situation
		const long double H = han+hbn+hcn+hdn -h(b+d,n) -h(a+c,n);
		if(H<0 && VERYCLOSE(H,0.0L))
			return 0.0L;
		assert(H>=0);
		return H;
	} else { // negatively-correlated, so we use the exception described in the paper.
		const long double H = h(c+d,n)+h(a+b,n);
		assert(H>=0);
		return H;
	}
}

long double H_X(in< vector<size_t> > X_sizes, const size_t n) {
	long double h_X = 0.0L;
	For(sz, *X_sizes) {
		h_X += h(*sz,n) + h(n-*sz,n);
	}
	return h_X;
}
long double calculate_H_X_given_Y( in< unordered_map< pair<size_t, size_t>, size_t > >intersections_, in< vector<size_t> > X_sizes, in< vector<size_t> > Y_sizes, const size_t N) {
	unordered_map< size_t, unordered_map<size_t, size_t> > intersections;
	For(inter, *intersections_) {
		intersections[inter->first.first][inter->first.second] = inter->second;
	}
	// Now, we can take each X_i, and try to find the Y_j which has the best H*(X_i|Y_j)
	long double H_X_Y = 0.0L;
	For(X_i, intersections) {
		const size_t i = X_i -> first;
		const size_t xi_size = X_sizes->at(i);
		in< unordered_map<size_t,size_t> > my_intersections = X_i->second;
		long double best_h_star_i_j_found_so_far = numeric_limits<long double> :: max();
		assert( best_h_star_i_j_found_so_far == numeric_limits<long double> :: max() );
		For(Y_j, *my_intersections) {
			const size_t j = Y_j -> first;
			const size_t yj_size = Y_sizes->at(j);

			const size_t d = Y_j -> second; // the intersection
			const size_t c = xi_size - d; // in X_i, but not in the intersection
			const size_t b = yj_size - d; // in Y_j, but not in the intersection
			const size_t a = N - d - c - b;
			assert(d>0);
			assert(xi_size >= d);
			assert(yj_size >= d);
			assert(N >= d+c+b);

			const long double h_star_i_j = H_star(a,b,c,d,N);
			assert(isfinite(h_star_i_j));
			assert(h_star_i_j >= 0);
			if(best_h_star_i_j_found_so_far > h_star_i_j)
				best_h_star_i_j_found_so_far = h_star_i_j;
		}
		// It's possible there were no intersections! In that case, define H(X_i|Y) = H(X_i)
		if(best_h_star_i_j_found_so_far == numeric_limits<long double> :: max()) {
			best_h_star_i_j_found_so_far = h(xi_size, N) + h(N-xi_size, N);
		}
		assert( best_h_star_i_j_found_so_far < numeric_limits<long double> :: max() );
		
		// OK, so that gets us H(X_i|Y), we need to add this into H(X|Y)
		H_X_Y += best_h_star_i_j_found_so_far;
	}
	return H_X_Y;
}
long double calculate_oNMI(lvalue_input :: in< std::vector< std::vector<int64_t> > > ground_truth, lvalue_input :: in<State> st) {
	// As described on arXiv http://arxiv.org/abs/1110.2515
	// I have to find the following:
	// 	H(X) and H(Y)
	// 	and I(X|Y), which will require the calculation of
	// 	H(X|Y) and H(Y|X)
	//
	// Naively, we would calculate H(X_i, Y_j) for every pair of communities
	// but we only need these for clusters where there is an overlap

	const int64_t N = st->N;

	if(st->get_K()==1) {
		return 0.0L; // We have just one cluster, and maybe some empty ones.  No NMI
	}

	// For each node, store (for both sides), the comms it is in
	vector< pair< vector<int>, vector<int> > > comms_this_node_is_in(N);
	vector<size_t> sizes_of_GT;
	vector<size_t> sizes_of_foundComms;
	for(size_t gt=0; gt<ground_truth->size(); ++gt) {
		sizes_of_GT.push_back( ground_truth->at(gt).size() );
		For(node, ground_truth->at(gt)) {
			comms_this_node_is_in.at(*node).first.push_back(gt);
		}
	}
	for(int k=0; k<st->get_K(); ++k) {
		vector<int64_t> my_nodes = st->get_comms().at(k).get_my_nodes_NO_COUNT();
		sizes_of_foundComms.push_back( my_nodes.size() );
		For(node, my_nodes) {
			comms_this_node_is_in.at(*node).second.push_back(k);
		}
	}
	//PP2(ground_truth->size(), st->get_K());
	assert(ground_truth->size() == sizes_of_GT.size());
	assert(st->get_K() == (int64_t)sizes_of_foundComms.size());

	unordered_map< pair<size_t, size_t>, size_t > intersections; // for each pair of comms with an intersection, record the size of it here
	for(int n=0; n<N; ++n) {
		For(left, comms_this_node_is_in.at(n).first) {
			For(right, comms_this_node_is_in.at(n).second) {
				intersections[ make_pair(*left, *right) ] ++;
			}
		}
	}


	// At this stage, intersections.first means ground truth, and intersections.second means found.
	// We'll swap them around later in the function

	long double H_GT = H_X(sizes_of_GT,N);
	long double H_GT_given_Fnd = calculate_H_X_given_Y(intersections, sizes_of_GT, sizes_of_foundComms, N);

	// Now, to swap the other way
	unordered_map< pair<size_t, size_t>, size_t > intersections_swapped;
	For(inter, intersections) {
		intersections_swapped[ make_pair(inter->first.second, inter->first.first) ] = inter->second;
	}

	long double H_Fnd = H_X(sizes_of_foundComms,N);
	if(VERYCLOSE(H_Fnd, 0.0L)) {
		return 0.0L; // We have just one cluster, and maybe some empty ones.  No NMI
	}
	long double H_Fnd_given_GT = calculate_H_X_given_Y(intersections_swapped, sizes_of_foundComms, sizes_of_GT, N);

	assert(H_GT  >= H_GT_given_Fnd || VERYCLOSE(H_GT , H_GT_given_Fnd));
	assert(H_Fnd >= H_Fnd_given_GT || VERYCLOSE(H_Fnd, H_Fnd_given_GT));

	const long double mutual_information = 0.5L * (H_GT + H_Fnd - H_GT_given_Fnd - H_Fnd_given_GT); // I_X_colon_Y
	const long double normalization_constant = H_GT > H_Fnd ? H_GT : H_Fnd;
	// assert(mutual_information <= normalization_constant);

	// PP(intersections.size());
	// For(inter, intersections) { cout << ' ' << inter->second; }
	// cout << endl;

	return mutual_information / normalization_constant;
}
