/*
 * - MK : add/remove empty cluster
 *
 * Partitioning of edges:
 * - Gibbs update of one edge
 * - SM-like move
 * - M3-like move
 *
 * Overlapping edges.
 * - Deprecated: M-H step on one edge.  Propose for each community independently, accept if and only if the row is non-empty
 * - Surely there is a Gibbs update I can do?
 * - SM-like edges - where there are *three* possible states, not just two
 */
#ifndef MOVES_HPP__
#define MOVES_HPP__

#include"score.hpp"
#include<utility>
#include <gsl/gsl_rng.h>
#include "lvalue_input.hpp"

void			seed_the_random_number_generator(int seed);
long double		gibbsUpdate(int64_t e, Score & sc);
long double		metroK(Score & sc);
std :: pair<long double,long double> 	gibbsUpdateJustTwoComms(
							int64_t e
							, Score & sc, const int64_t main_cluster
							, const int64_t secondary_cluster
							, std :: pair<int,int> possibly_force = std :: make_pair(-1,-1)
							);
long double		split_or_merge(Score & sc);
long double		M3(Score & sc);
long double		split_or_merge_by_seed_expansion(Score & sc);
long double		split_or_merge_on_a_shared_edge(Score & sc);
long double one_node_simple_update(Score &sc);
long double		one_node_SIMPLEST_update(Score &sc, int64_t n);
long double		gibbs_one_comm_one_edge(Score & sc, const int64_t e);
long double 		gibbsUpdateNearby(Score & sc, int64_t e);
void                    expand_seed(const int comm_source, const int comm_new, const int seed_edge, const std:: vector<int> &E, Net net, long double &delta_score, Score &sc);
std:: vector<int>             my_edges_in_a_random_order(State &st, const int k);

gsl_rng * rng();

// Some constraints enforced from outside, only for the first 100 iterations or so
extern int64_t GLOBAL_constraint_min_K;

#endif
