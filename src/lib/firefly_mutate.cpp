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
#include "kmeans.h"

sz firefly_count_mutable_entities(const conf &c)
{
    sz counter = 0;
    VINA_FOR_IN(i, c.ligands)
    counter += 2 + c.ligands[i].torsions.size();
    VINA_FOR_IN(i, c.flex)
    counter += c.flex[i].torsions.size();
    return counter;
}

void update_2_and_3(int i, int order0, firefly *fireflies, output_type &candidate_2, output_type &candidate_3, output_type &candidate) 
{
    KMeans kmeans(3, 100);
    double  all_points[fireflies->number];
    int * cluster_best = new int[3];
    for (int k = 0; k < fireflies->number; k++)
    {
        all_points[k] = fireflies->getCurrentFit(k);
    }
        
    kmeans.run(all_points, fireflies->number);
    cluster_best[0] = kmeans.a;
    cluster_best[1] = kmeans.b;
    cluster_best[2] = kmeans.c;
    
    int count = 0;
    int other_cluster[2];
    for (int k = 0; k < 3; k++ )
    {
        if (cluster_best[k]  != order0)
        {    
            other_cluster[count] = cluster_best[k];
            count++;
        } 
    }
    int second, third;

    if (fireflies->getCurrentFit(other_cluster[0]) < fireflies->getCurrentFit(other_cluster[1]))
    {
        second = other_cluster[0];
        third = other_cluster[1];
    }
    else
    {
        second = other_cluster[1];
        third = other_cluster[0];
    }

    for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
        candidate_2.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(second, z);
    candidate_2.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(second);
    candidate_2.c.ligands[i].rigid.position = fireflies->getCurrentPosition(second);
    candidate_2.e = fireflies->getCurrentFit(second);

    for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
        candidate_3.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(third, z);
    candidate_3.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(third);
    candidate_3.c.ligands[i].rigid.position = fireflies->getCurrentPosition(third);
    candidate_3.e = fireflies->getCurrentFit(third);
}

void firefly_mutate_conf(
    output_type &candidate, 
    output_type &candidate_1, 
    output_type &candidate_2, 
    output_type &candidate_3, 
    const model &m, 
    fl amplitude, 
    rng &generator, 
    firefly *fireflies, 
    double *PersonalBest, 
    const precalculate &p, 
    const igrid &ig, 
    change &g, 
    const vec &v, 
    quasi_newton &quasi_newton_par, 
    int step
    )
{
    int shrink = 5; //roughing factor R = shrink / 10
    int cr = 20;    //roughing condition
    output_type tmp_2 = candidate;
    //if (step == 0)
        //std::cout << "current_alpha:" << fireflies->alpha << '\n';
    
    //chaotic 
    double temp = fireflies->gamma;
    if (fireflies->chaos == 1)  
        fireflies->gamma = fireflies->mu1 * fireflies->gamma * (1-fireflies->gamma);
    
    if (fireflies->chaos == 1 && fireflies->levy_flight == 0)
        fireflies->alpha = fireflies->mu2 * temp * (1-fireflies->alpha);

    //fireflies->beta = fireflies->mu1 * fireflies->beta * (1-fireflies->beta);
    //fireflies->beta = std::exp(1) - 4.9 * fireflies->beta * fireflies->beta - 0.58;
    //fireflies->gamma = std::exp(1) - 4.9 * fireflies->gamma * fireflies->gamma - 0.58;
    //fireflies->alpha = fireflies->mu2 * temp * (1-fireflies->alpha);
    //std::cout << "current_step:" << step << '\n';
    fireflies->alpha = fireflies->alpha * 0.99;
    //fireflies->alpha = 0.9 - 0.9 * pow((step-1)/1999, 2.1);
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
                    candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(y, z);

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
                        fireflies->updateGlobalBestFit(tmp_2.e);
                        fireflies->updateGlobalBestFirefly(y);
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

            if (fireflies->clustering == 1)
                update_2_and_3(i, order[0], fireflies, candidate_2, candidate_3, candidate);

            for (int k = 0; k < fireflies->number; k++)
                for (int j = 0; j < fireflies->number; j++)
                    if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                    {
                        fireflies->moveFireflyPosition(order[j], k, generator);
                    }

            
            if (fireflies->elite == 0)
                fireflies->moveFireflyPositionRandomly(order[0], generator);
            else
            {
                vec original_position = fireflies->getCurrentPosition(order[0]);
                bool success_pos = false;

                for (int x = 0; x < 20; x++)
                {
                    fireflies->updateCurrentPosition(order[0], original_position);
                    fireflies->moveFireflyPositionRandomly(order[0], generator);

                    output_type tmp_3 = candidate;

                    tmp_3.c.ligands[i].rigid.position = fireflies->getCurrentPosition(order[0]);
                    tmp_3.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(order[0]);
                    for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                        tmp_3.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(order[0], z);

                    quasi_newton_par(tmp_m, p, ig, tmp_3, g, v, true);

                    if (tmp_3.e < firefly::gbest_fit)
                    {
                        PersonalBest[order[0]] = tmp_3.e;
                        fireflies->updateCurrentFit(order[0], tmp_3.e);
                        fireflies->updateCurrentPosition(order[0], tmp_3.c.ligands[i].rigid.position);
                        fireflies->updateCurrentOrientation(order[0], tmp_3.c.ligands[i].rigid.orientation);
                        for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                            fireflies->updateCurrentTorsion(order[0], tmp_3.c.ligands[i].torsions[z], z);

                        fireflies->updateGlobalBestFit(tmp_3.e);

                        firefly::gbest_position = tmp_3.c.ligands[i].rigid.position;
                        firefly::gbest_orientation = tmp_3.c.ligands[i].rigid.orientation;
                        for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                            firefly::gbest_torsion[z] = tmp_3.c.ligands[i].torsions[z];
                        
                        success_pos = true;
                        break;
                    }
                }

                if (success_pos == false)
                {
                    fireflies->updateCurrentPosition(order[0], original_position);
                }
            }

            

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
                        candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(y, z);

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
                            fireflies->updateCurrentFit(y, tmp_2.e);
                            fireflies->updateCurrentPosition(y, tmp_2.c.ligands[i].rigid.position);
                            fireflies->updateCurrentOrientation(y, tmp_2.c.ligands[i].rigid.orientation);
                            for (int z = 0; z < tmp_2.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(y, tmp_2.c.ligands[i].torsions[z], z);

                            //set the global best(energy value and position);
                           
                        }

                        if (tmp_2.e < firefly::gbest_fit)
                        {

                            fireflies->updateGlobalBestFit(tmp_2.e);
                            fireflies->updateGlobalBestFirefly(y);
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

                if (fireflies->clustering == 1)
                    update_2_and_3(i, order[0], fireflies, candidate_2, candidate_3, candidate);\

                for (int k = 0; k < fireflies->number; k++)
                    for (int j = 0; j < fireflies->number; j++)
                        if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                        {
                            fireflies->moveFireflyOrientation(order[j], k, generator);
                        }

                if (fireflies->elite == 0)
                    fireflies->moveFireflyOrientationRandomly(order[0], generator);
                else
                {
                    qt original_orientation = fireflies->getCurrentOrientation(order[0]);
                    bool success_ori = false;

                    for (int x = 0; x < 20; x++)
                    {
                        fireflies->updateCurrentOrientation(order[0], original_orientation);
                        fireflies->moveFireflyOrientationRandomly(order[0], generator);

                        output_type tmp_3 = candidate;

                        tmp_3.c.ligands[i].rigid.position = fireflies->getCurrentPosition(order[0]);
                        tmp_3.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(order[0]);
                        for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                            tmp_3.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(order[0], z);

                        quasi_newton_par(tmp_m, p, ig, tmp_3, g, v, true);

                        if (tmp_3.e < firefly::gbest_fit)
                        {
                            PersonalBest[order[0]] = tmp_3.e;
                            fireflies->updateCurrentFit(order[0], tmp_3.e);
                            fireflies->updateCurrentPosition(order[0], tmp_3.c.ligands[i].rigid.position);
                            fireflies->updateCurrentOrientation(order[0], tmp_3.c.ligands[i].rigid.orientation);
                            for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(order[0], tmp_3.c.ligands[i].torsions[z], z);

                            fireflies->updateGlobalBestFit(tmp_3.e);

                            firefly::gbest_position = tmp_3.c.ligands[i].rigid.position;
                            firefly::gbest_orientation = tmp_3.c.ligands[i].rigid.orientation;
                            for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                                firefly::gbest_torsion[z] = tmp_3.c.ligands[i].torsions[z];
                            
                            success_ori = true;
                            break;
                        }
                    }

                    if (success_ori == false)
                    {
                        fireflies->updateCurrentOrientation(order[0], original_orientation);
                    }
                }

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
                        candidate.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(y, z);

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
                            for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(y, tmp_2.c.ligands[i].torsions[z], z);

                            //set the global best(energy value and position);
                            
                        }

                        if (tmp_2.e < firefly::gbest_fit)
                        {
                            fireflies->updateGlobalBestFit(tmp_2.e);
                            fireflies->updateGlobalBestFirefly(y);
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

                if (fireflies->clustering == 1)
                    update_2_and_3(i, order[0], fireflies, candidate_2, candidate_3, candidate);

                for (int k = 0; k < fireflies->number; k++)
                    for (int j = 0; j < fireflies->number; j++)
                        if (fireflies->getCurrentFit(k) > fireflies->getCurrentFit(order[j]))
                            fireflies->moveFireflyTorsion(order[j], k, generator, which);

                if (fireflies->elite == 0)
                    fireflies->moveFireflyTorsionRandomly(order[0], generator, which);
                else
                {
                    fl original_torsion = fireflies->getCurrentTorsion(order[0], which);
                    bool success_tor = false;
                    
                    for (int x = 0; x < 20; x++)
                    {
                        fireflies->updateCurrentTorsion(order[0], original_torsion, which);
                        fireflies->moveFireflyTorsionRandomly(order[0], generator, which);

                        output_type tmp_3 = candidate;

                        tmp_3.c.ligands[i].rigid.position = fireflies->getCurrentPosition(order[0]);
                        tmp_3.c.ligands[i].rigid.orientation = fireflies->getCurrentOrientation(order[0]);
                        for (int z = 0; z < candidate.c.ligands[i].torsions.size(); z++)
                            tmp_3.c.ligands[i].torsions[z] = fireflies->getCurrentTorsion(order[0], z);

                        quasi_newton_par(tmp_m, p, ig, tmp_3, g, v, true);

                        if (tmp_3.e < firefly::gbest_fit)
                        {
                            PersonalBest[order[0]] = tmp_3.e;
                            fireflies->updateCurrentFit(order[0], tmp_3.e);
                            fireflies->updateCurrentPosition(order[0], tmp_3.c.ligands[i].rigid.position);
                            fireflies->updateCurrentOrientation(order[0], tmp_3.c.ligands[i].rigid.orientation);
                            for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                                fireflies->updateCurrentTorsion(order[0], tmp_3.c.ligands[i].torsions[z], z);

                            fireflies->updateGlobalBestFit(tmp_3.e);

                            firefly::gbest_position = tmp_3.c.ligands[i].rigid.position;
                            firefly::gbest_orientation = tmp_3.c.ligands[i].rigid.orientation;
                            for (int z = 0; z < tmp_3.c.ligands[i].torsions.size(); z++)
                                firefly::gbest_torsion[z] = tmp_3.c.ligands[i].torsions[z];
                            
                            success_tor = true;
                            break;
                        }
                    }

                    if (success_tor == false)
                    {
                        fireflies->updateCurrentTorsion(order[0], original_torsion, which);
                    }
                }
                
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
