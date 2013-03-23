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
std :: pair<long double,long double> 	gibbsUpdateJustTwoComms(int64_t e, Score & sc, const int64_t main_cluster, const int64_t secondary_cluster) ;

#endif
