#ifndef VINA_PARALLEL_FIREFLY_H
#define VINA_PARALLEL_FIREFLY_H

#include "firefly_search.h"

struct parallel_firefly
{
        firefly_search firefly;
        sz num_tasks;
        sz num_threads;
        bool display_progress;
        parallel_firefly() : num_tasks(8), num_threads(1), display_progress(true) {}
        void operator()(const model &m,
                        output_container &out,
                        output_container &out_2,
                        output_container &out_3,
                        const precalculate &p,
                        const igrid &ig,
                        const precalculate &p_widened,
                        const igrid &ig_widened,
                        const vec &corner1,
                        const vec &corner2,
                        rng &generatorint,
                        int num_of_fireflies,
                        double gamma,
                        double beta,
                        double alpha,
                        double mu1, 
                        double mu2) const;
};

#endif
