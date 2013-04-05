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
long double		split_or_merge_on_a_shared_edge(Score & sc);
long double one_node_simple_update(Score &sc);
long double		one_node_SIMPLEST_update(Score &sc, int64_t n);

#endif
