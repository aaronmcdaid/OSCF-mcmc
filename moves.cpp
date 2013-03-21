#include "moves.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cassert>

gsl_rng * r = NULL;
void			seed_the_random_number_generator(int seed) {
					assert( r == NULL );
					r = gsl_rng_alloc (gsl_rng_taus);
					gsl_rng_set(r, seed);
}
