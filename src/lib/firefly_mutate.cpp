/*
        FireflyVina version 1.0                     

        This file is revised from mutate.cpp in AutoDock Vina.

        Authors: Marcus M. C. Ao  <mingchi_@hotmail.com>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "firefly_mutate.h"

sz firefly_count_mutable_entities(const conf &c)
{
    sz counter = 0;
    VINA_FOR_IN(i, c.ligands)
    counter += 2 + c.ligands[i].torsions.size();
    VINA_FOR_IN(i, c.flex)
    counter += c.flex[i].torsions.size();
    return counter;
}

void firefly_mutate_conf(output_type &candidate, output_type &candidate_1, const model &m, fl amplitude, rng &generator, firefly *fireflies, double *PersonalBest, const precalculate &p, const igrid &ig, change &g, const vec &v, quasi_newton &quasi_newton_par, int step)
{
    int shrink = 5; //roughing factor R = shrink / 10
    int cr = 20;    //roughing condition
    output_type tmp_2 = candidate;
    //if (step == 0)
        //std::cout << "current_alpha:" << fireflies->alpha << '\n';
    
    //chaotic 
    double temp = fireflies->gamma;   
    fireflies->gamma = fireflies->mu1 * fireflies->gamma * (1-fireflies->gamma);
    
    fireflies->alpha = fireflies->mu2 * temp * (1-fireflies->alpha);
    
    //fireflies->alpha = fireflies->alpha * 0.97;
    //std::cout << "current_beta:" << fireflies->beta << '\n';
    //std::cout << "current_alpha:" << fireflies->alpha << '\n';
    //std::cout << "current_number:" << fireflies->number << '\n';
    sz mutable_entities_num = firefly_count_mutable_entities(candidate.c);
    if (mutable_entities_num == 0)
        return;
    int which_int = random_int(0, int(mutable_entities_num - 1), generator);
    VINA_CHECK(which_int >= 0);
    sz which = sz(which_int);
    VINA_CHECK(which < mutable_entities_num);

    int y;
    VINA_FOR_IN(i, candidate.c.ligands)
    {
        model tmp_m = m;
        const vec authentic_v(1000, 1000, 1000);

        //Take part the position (either position or orientation or torsion)
        if (which == 0)
        {

            //loop for each firefly
            for (y = 0; y < fireflies->number; y++)
            {

                candidate.c.ligands[i].rigid.position = fireflies->getCurrentPosition(y);
                candidate.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(y);
                for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                    candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(i, z);

                //local search
                quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

                //set the personal best(energy value and position);
                if (candidate.e < fireflies->getPersonalBest(y) || step <= cr)
                {
                    if (candidate.e < fireflies->getPersonalBest(y))
                        fireflies->updatePersonalBest(y, candidate.e);
                    tmp_2 = candidate;

                    quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);

                    if (tmp_2.e < PersonalBest[y])
                    {
                        PersonalBest[y] = tmp_2.e;

                        fireflies->updateCurrentFit(y, tmp_2.e);

                        fireflies->updateCurrentPosition(y, tmp_2.c.ligands[i].rigid.position);
                        fireflies->updateCurrentOrientation(y, tmp_2.c.ligands[i].rigid.orientation);
                        for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                            fireflies->updateCurrentTorsion(y, tmp_2.c.ligands[i].torsions[z], z);
                    }

                    //set the global best(energy value and position);
                    if (tmp_2.e < firefly::gbest_fit)
                    {
                        //std::cout << "current_P:" << candidate.e << "	quasi_e:"<<candidate.e << "	the current best_P:" << firefly::gbest_fit<<"	Current number steps" <<step <<'\n';
                        fireflies->updateGlobalBestFit(tmp_2.e);

                        firefly::gbest_position = tmp_2.c.ligands[i].rigid.position;
                        firefly::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
                        for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                            firefly::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
                    }
                }
            }

            int order[fireflies->number];
            for (int x = 0; x < fireflies->number; x++)
                order[x] = x;

            for (int t = 0; t < fireflies->number - 1; t++)
                for (int u = t + 1; u < fireflies->number; u++)
                    if (fireflies->getCurrentFit(order[t]) > fireflies->getCurrentFit(order[u]))
                    {
                        int temp = order[t];
                        order[t] = order[u];
                        order[u] = temp;
                    }

            for (int k = 0; k < fireflies->number; k++)
                for (int j = 0; j < fireflies->number; j++)
                    if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                    {
                        fireflies->moveFireflyPosition(order[j], k, generator);
                    }

            fireflies->moveFireflyPositionRandomly(order[0], generator);

            for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                candidate_1.c.ligands[i].torsions[z] = firefly::gbest_torsion[z];
            candidate_1.c.ligands[i].rigid.orientation = firefly::gbest_orientation;
            candidate_1.c.ligands[i].rigid.position = firefly::gbest_position;
            candidate_1.e = firefly::gbest_fit;
            return;
        }

        --which;
        //Take part orientation (either position or orientation or torsion)
        if (which == 0)
        {
            fl gr = m.gyration_radius(i);
            if (gr > epsilon_fl)
            { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?

                for (y = 0; y < fireflies->number; y++)
                {

                    candidate.c.ligands[i].rigid.position = fireflies->getCurrentPosition(y);
                    candidate.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(y);
                    for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                        candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(i, z);

                    quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

                    //set the personal best(energy value and position);
                    if (candidate.e < fireflies->getPersonalBest(y) || step <= cr)
                    {

                        if (candidate.e < fireflies->getPersonalBest(y))
                            fireflies->updatePersonalBest(y, candidate.e);

                        tmp_2 = candidate;

                        //full local search
                        quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);

                        if (tmp_2.e < PersonalBest[y])
                        {
                            PersonalBest[y] = tmp_2.e;
                            fireflies->updateCurrentPosition(y, tmp_2.c.ligands[i].rigid.position);
                            fireflies->updateCurrentOrientation(y, tmp_2.c.ligands[i].rigid.orientation);
                            for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(y, tmp_2.c.ligands[i].torsions[z], z);

                            //set the global best(energy value and position);
                            if (tmp_2.e < firefly::gbest_fit)
                            {

                                fireflies->updateGlobalBestFit(tmp_2.e);

                                firefly::gbest_position = tmp_2.c.ligands[i].rigid.position;
                                firefly::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
                                for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                                    firefly::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
                            }
                        }
                    }
                }

                int order[fireflies->number];
                for (int x = 0; x < fireflies->number; x++)
                    order[x] = x;

                for (int t = 0; t < fireflies->number - 1; t++)
                    for (int u = t + 1; u < fireflies->number; u++)
                        if (fireflies->getCurrentFit(order[t]) > fireflies->getCurrentFit(order[u]))
                        {
                            int temp = order[t];
                            order[t] = order[u];
                            order[u] = temp;
                        }

                for (int k = 0; k < fireflies->number; k++)
                    for (int j = 0; j < fireflies->number; j++)
                        if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                        {
                            fireflies->moveFireflyOrientation(order[j], k, generator);
                        }
                
                fireflies->moveFireflyOrientationRandomly(order[0], generator);

                for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                    candidate_1.c.ligands[i].torsions[z] = firefly::gbest_torsion[z];
                candidate_1.c.ligands[i].rigid.orientation = firefly::gbest_orientation;
                candidate_1.c.ligands[i].rigid.position = firefly::gbest_position;
                candidate_1.e = firefly::gbest_fit;
                return;
            }

            /*Torsions*/
            --which;
            if (which < candidate.c.ligands[i].torsions.size())
            {

                for (y = 0; y < fireflies->number; y++)
                {

                    candidate.c.ligands[i].rigid.position = fireflies->getCurrentPosition(y);
                    candidate.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(y);
                    for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                        candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(i, z);

                    quasi_newton_par(tmp_m, p, ig, candidate, g, v, shrink);

                    //set the personal best(energy value and position);
                    if (candidate.e < fireflies->getPersonalBest(y) || step <= cr)
                    {
                        if (candidate.e < fireflies->getPersonalBest(y))
                            fireflies->updatePersonalBest(y, candidate.e);

                        tmp_2 = candidate;
                        quasi_newton_par(tmp_m, p, ig, tmp_2, g, v);

                        if (tmp_2.e < PersonalBest[y])
                        {
                            fireflies->updateCurrentPosition(y, tmp_2.c.ligands[i].rigid.position);
                            fireflies->updateCurrentOrientation(y, tmp_2.c.ligands[i].rigid.orientation);
                            for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(y, tmp_2.c.ligands[i].torsions[z], z);

                            //set the global best(energy value and position);
                            if (tmp_2.e < firefly::gbest_fit)
                            {
                                //std::cout << "current_T:" << candidate.e << "	quasi_e:"<<candidate.e<< "	the current best_T:" << firefly::gbest_fit<<"	Current number step"<<step <<'\n';
                                fireflies->updateGlobalBestFit(tmp_2.e);

                                firefly::gbest_position = tmp_2.c.ligands[i].rigid.position;
                                firefly::gbest_orientation = tmp_2.c.ligands[i].rigid.orientation;
                                for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                                    firefly::gbest_torsion[z] = tmp_2.c.ligands[i].torsions[z];
                            }
                        }
                    }
                }

                int order[fireflies->number];
                for (int x = 0; x < fireflies->number; x++)
                    order[x] = x;

                for (int t = 0; t < fireflies->number - 1; t++)
                    for (int u = t + 1; u < fireflies->number; u++)
                        if (fireflies->getCurrentFit(order[t]) > fireflies->getCurrentFit(order[u]))
                        {
                            int temp = order[t];
                            order[t] = order[u];
                            order[u] = temp;
                        }

                for (int k = 0; k < fireflies->number; k++)
                    for (int j = 0; j < fireflies->number; j++)
                        if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                        {
                            fireflies->moveFireflyTorsion(order[j], k, generator, which);
                        }

                fireflies->moveFireflyTorsionRandomly(order[0], generator, which);

                for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                    candidate_1.c.ligands[i].torsions[z] = firefly::gbest_torsion[z];
                candidate_1.c.ligands[i].rigid.orientation = firefly::gbest_orientation;
                candidate_1.c.ligands[i].rigid.position = firefly::gbest_position;
                candidate_1.e = firefly::gbest_fit;
                return;
            }
            which -= candidate.c.ligands[i].torsions.size();
        }

        VINA_FOR_IN(i, candidate.c.flex)
        {
            if (which < candidate.c.flex[i].torsions.size())
            {
                candidate.c.flex[i].torsions[which] = random_fl(-pi, pi, generator);
                return;
            }
            which -= candidate.c.flex[i].torsions.size();
        }
    }
}
