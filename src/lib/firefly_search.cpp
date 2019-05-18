#include "firefly_search.h"
#include "coords.h"
#include "quasi_newton.h"
#include "firefly_mutate.h"
#include "mutate.h"

output_type firefly_search::operator()(model &m,
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
                                       double lambda) const
{
    output_container tmp;
    this->operator()(m,
                     tmp,
                     p,
                     ig,
                     p_widened,
                     ig_widened,
                     corner1,
                     corner2,
                     increment_me,
                     generator,
                     num_of_fireflies,
                     gamma,
                     beta,
                     alpha,
                     mu1,
                     mu2,
                     lambda); // call the version that produces the whole container
    VINA_CHECK(!tmp.empty());
    return tmp.front();
}

bool metropolis_accept(fl old_f, fl new_f, fl temperature, rng &generator)
{
    if (new_f < old_f)
        return true;
    const fl acceptance_probability = std::exp((old_f - new_f) / temperature);
    return random_fl(0, 1, generator) < acceptance_probability;
}

// out is sorted
void firefly_search::operator()(model &m,
                                output_container &out,
                                const precalculate &p,
                                const igrid &ig,
                                const precalculate &p_widened,
                                const igrid &ig_widened,
                                const vec &corner1,
                                const vec &corner2,
                                incrementable *increment_me,
                                rng &generator,
                                int num_fireflies,
                                double gamma,
                                double beta,
                                double alpha,
                                double mu1,
                                double mu2,
                                double lambda) const
{
    vec authentic_v(1000, 1000, 1000); // FIXME? this is here to avoid max_fl/max_fl
    conf_size s = m.get_size();
    change g(s);
    output_type tmp(s, 0);
    tmp.c.randomize(corner1, corner2, generator); //first randomize
    fl best_e = max_fl;
    quasi_newton quasi_newton_par;
    quasi_newton_par.max_steps = ssd_par.evals;
    output_type tmp_rough = tmp;
    firefly fireflies(num_fireflies, gamma, beta, alpha, mu1, mu2, lambda, corner1, corner2, generator, tmp.c);
    double* PersonalBest = new double[1000];
    for(int cou = 0; cou < 100; cou++)
		PersonalBest[cou] = 0;
    double energy = 0;
    int count = 0;

    //printf("MAXSTEP %d\n", num_steps);
    //printf("TORSIZE %d\n", tmp.c.ligands[0].torsions.size());
    //printf("PARTICLE %d\n", particle.number);

    VINA_U_FOR(step, num_steps)
    {

        if (increment_me)
            ++(*increment_me);
        output_type candidate = tmp_rough;
        output_type candidate_1 = tmp;

        firefly_mutate_conf(candidate, candidate_1, m, mutation_amplitude, generator, &fireflies, PersonalBest, p, ig, g, hunt_cap, quasi_newton_par, step); //for each particle loop

        tmp_rough = candidate;
        if (step == 0 || metropolis_accept(tmp.e, candidate_1.e, temperature, generator))
        {
            tmp = candidate_1;
            m.set(tmp.c); // FIXME? useless?

            // FIXME only for very promising ones
            if (tmp.e < best_e || out.size() < num_saved_mins)
            {
                quasi_newton_par(m, p, ig, tmp, g, authentic_v);
                m.set(tmp.c); // FIXME? useless?
                tmp.coords = m.get_heavy_atom_movable_coords();
                add_to_output_container(out, tmp, min_rmsd, num_saved_mins); // 20 - max size
                if (tmp.e < best_e)
                    best_e = tmp.e;
            }
        }

        /***Criteria defined by PSOVina***/

        if (std::abs(firefly::gbest_fit - energy) < 0.0001)
        {

            count += 1;
            if (count > 350)
            {
                //printf("Terminated: %d \n",step);
                //printf("CONVERGEAT %d\n", step);
                step = num_steps; //break the loop
                count = 0;
            }
        }
        else
        {
            energy = firefly::gbest_fit;
            count = 0;
        }
    }

    //printf("GBEST %lf\n", firefly::gbest_fit);

    VINA_CHECK(!out.empty());
    VINA_CHECK(out.front().e <= out.back().e); // make sure the sorting worked in the correct order
}
