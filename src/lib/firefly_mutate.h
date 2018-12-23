#ifndef VINA_FIREFLY_MUTATE_H
#define VINA_FIREFLY_MUTATE_H

#include "firefly.h"
#include "model.h"
#include "quasi_newton.h"

// does not set model
void firefly_mutate_conf(output_type &c,
                         output_type &c_1,
                         const model &m,
                         fl amplitude,
                         rng &generator,
                         firefly *,
                         double *,
                         const precalculate &,
                         const igrid &,
                         change &,
                         const vec &,
                         quasi_newton &,
                         int);

void firefly_mutate_conf(output_type &c,
                         const model &m,
                         fl amplitude,
                         rng &generator,
                         firefly *,
                         const precalculate &,
                         const igrid &,
                         change &,
                         const vec &,
                         quasi_newton &,
                         int);
#endif
