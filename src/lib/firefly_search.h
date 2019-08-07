#ifndef VINA_MONTE_CARLO_H
#define VINA_MONTE_CARLO_H

#include "ssd.h"
#include "incrementable.h"
#include "firefly.h"

struct firefly_search
{
    unsigned num_steps;
    fl temperature;
    vec hunt_cap;
    fl min_rmsd;
    sz num_saved_mins;
    fl mutation_amplitude;
    ssd ssd_par;

    firefly_search() : num_steps(2500), temperature(1.2), hunt_cap(10, 1.5, 10), min_rmsd(0.5), num_saved_mins(50), mutation_amplitude(2) {} // T = 600K, R = 2cal/(K*mol) -> temperature = RT = 1.2;  num_steps = 50*lig_atoms = 2500

    output_type operator()(model &m,
                           const precalculate &p,
                           const igrid &ig,
                           const precalculate &p_widened,
                           const igrid &ig_widened,
                           const vec &corner1,
                           const vec &corner2,
                           incrementable *increment_me,
                           rng &generator,
                           int num_of_fireflies,
                           double gamma,
                           double beta,
                           double alpha,
                           double mu1,
                           double mu2,
                           double lambda,
                           int clustering,
                           int levy_flight,
                           int chaos,
                           int elite) const;
    // out is sorted
    void operator()(model &m,
                    output_container &out,
                    output_container &out_2,
                    output_container &out_3,
                    const precalculate &p,
                    const igrid &ig,
                    const precalculate &p_widened,
                    const igrid &ig_widened,
                    const vec &corner1,
                    const vec &corner2,
                    incrementable *increment_me,
                    rng &generator,
                    int num_of_fireflies,
                    double gamma,
                    double beta,
                    double alpha,
                    double mu1,
                    double mu2,
                    double lambda,
                    int clustering,
                    int levy_flight,
                    int chaos,
                    int elite) const;
};

#endif
