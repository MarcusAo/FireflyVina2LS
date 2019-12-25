#include "parallel.h"
#include "parallel_firefly.h"
#include "coords.h"
#include "parallel_progress.h"

struct parallel_firefly_task
{
    model m;
    output_container out;
    output_container out_2;
    output_container out_3;
    rng generator;
    parallel_firefly_task(const model &m_, int seed) : m(m_), generator(static_cast<rng::result_type>(seed)) {}
};

typedef boost::ptr_vector<parallel_firefly_task> parallel_firefly_task_container;

struct parallel_firefly_aux
{
    const firefly_search *firefly;
    const precalculate *p;
    const igrid *ig;
    const precalculate *p_widened;
    const igrid *ig_widened;
    const vec *corner1;
    const vec *corner2;
    const int num_of_fireflies, clustering, levy_flight, chaos, elite;
    const double gamma, beta, alpha, mu1, mu2, lambda;
    parallel_progress *pg;
    parallel_firefly_aux(const firefly_search *firefly_,
                         const precalculate *p_,
                         const igrid *ig_,
                         const precalculate *p_widened_,
                         const igrid *ig_widened_,
                         const vec *corner1_,
                         const vec *corner2_,
                         parallel_progress *pg_,
                         int num_of_fireflies_,
                         double gamma_,
                         double beta_,
                         double alpha_,
                         double mu1_,
                         double mu2_,
                         double lambda_,
                         int clustering_,
                         int levy_flight_,
                         int chaos_,
                         int elite_) : firefly(firefly_),
                                       p(p_),
                                       ig(ig_),
                                       p_widened(p_widened_),
                                       ig_widened(ig_widened_),
                                       corner1(corner1_),
                                       corner2(corner2_), pg(pg_),
                                       num_of_fireflies(num_of_fireflies_),
                                       gamma(gamma_),
                                       beta(beta_),
                                       alpha(alpha_),
                                       mu1(mu1_),
                                       mu2(mu2_),
                                       lambda(lambda_),
                                       clustering(clustering_),
                                       levy_flight(levy_flight_),
                                       chaos(chaos_),
                                       elite(elite_) {}

    void operator()(parallel_firefly_task &t) const
    {
        (*firefly)(t.m,
                   t.out,
                   t.out_2,
                   t.out_3,
                   *p,
                   *ig,
                   *p_widened,
                   *ig_widened,
                   *corner1,
                   *corner2,
                   pg,
                   t.generator,
                   num_of_fireflies,
                   gamma,
                   beta,
                   alpha,
                   mu1,
                   mu2,
                   lambda,
                   clustering,
                   levy_flight,
                   chaos,
                   elite);
    }
};

void merge_output_containers(const output_container &in, output_container &out, fl min_rmsd, sz max_size)
{
    VINA_FOR_IN(i, in)
        add_to_output_container(out, in[i], min_rmsd, max_size);
}

void merge_output_containers(const parallel_firefly_task_container &many, output_container &out, fl min_rmsd, sz max_size)
{
    min_rmsd = 0.000001; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    VINA_FOR_IN(i, many)
    merge_output_containers(many[i].out, out, min_rmsd, max_size);
    out.sort();
}

void merge_output_containers2(const parallel_firefly_task_container &many, output_container &out, fl min_rmsd, sz max_size)
{
    min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    VINA_FOR_IN(i, many)
    merge_output_containers(many[i].out_2, out, min_rmsd, max_size);
    out.sort();
}

void merge_output_containers3(const parallel_firefly_task_container &many, output_container &out, fl min_rmsd, sz max_size)
{
    min_rmsd = 2; // FIXME? perhaps it's necessary to separate min_rmsd during search and during output?
    VINA_FOR_IN(i, many)
    merge_output_containers(many[i].out_3, out, min_rmsd, max_size);
    out.sort();
}

void parallel_firefly::operator()(const model &m,
                                  output_container &out,
                                  output_container &out_2,
                                  output_container &out_3,
                                  const precalculate &p,
                                  const igrid &ig,
                                  const precalculate &p_widened,
                                  const igrid &ig_widened,
                                  const vec &corner1,
                                  const vec &corner2,
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
                                  int elite) const
{
    parallel_progress pp;
    parallel_firefly_aux parallel_firefly_aux_instance(
        &firefly,
        &p,
        &ig,
        &p_widened,
        &ig_widened,
        &corner1,
        &corner2,
        (display_progress ? (&pp) : NULL),
        num_of_fireflies,
        gamma,
        beta,
        alpha,
        mu1,
        mu2,
        lambda,
        clustering,
        levy_flight,
        chaos,
        elite);

    parallel_firefly_task_container task_container;
    VINA_FOR(i, num_tasks)
    task_container.push_back(new parallel_firefly_task(m, random_int(0, 1000000, generator)));
    if (display_progress)
        pp.init(num_tasks * firefly.num_steps);
    parallel_iter<parallel_firefly_aux, parallel_firefly_task_container, parallel_firefly_task, true> parallel_iter_instance(&parallel_firefly_aux_instance, num_threads);
    parallel_iter_instance.run(task_container);
    merge_output_containers(task_container, out, firefly.min_rmsd, firefly.num_saved_mins);
    if (clustering == 1)
    {
        merge_output_containers2(task_container, out_2, firefly.min_rmsd, firefly.num_saved_mins);
        merge_output_containers3(task_container, out_3, firefly.min_rmsd, firefly.num_saved_mins);
    }
}
