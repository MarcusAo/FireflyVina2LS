#include "firefly.h"
#include "random.h"

// the vector of degree of freedom
qt firefly::gbest_orientation;
fl *firefly::gbest_torsion;
vec firefly::gbest_position;

double firefly::gbest_fit;

firefly::firefly(int num_fireflies, double gamma, double beta, double alpha, const vec corner1, const vec corner2, rng &g, conf &c)
{
    sz torsionSize = c.ligands[0].torsions.size();
    this->gamma = gamma,
    this->beta = beta,
    this->alpha = alpha;
    this->number = num_fireflies;
    this->g = g;
    this->corner1[0] = corner1[0]; //minmum
    this->corner1[1] = corner1[1];
    this->corner1[2] = corner1[2];
    this->corner2[0] = corner2[0]; //maximum
    this->corner2[1] = corner2[1];
    this->corner2[2] = corner2[2];
    this->torsionSize = (int)torsionSize;
    this->R1Max_ = 1;
    this->R1Min_ = 0;
    this->R2Max_ = 1;
    this->R2Min_ = 0;
    firefly::gbest_torsion = new fl[torsionSize];
    init(g, c);
}

void firefly::init(rng &g, conf &c)
{
    int i;
    for (i = 0; i < this->number; i++)
    {
        lampyridae single_lampyridae;

        single_lampyridae.current_fit = 1.7976931348623158e+308;
        single_lampyridae.pbest_fit = 1.7976931348623158e+308;

        //set position part
        single_lampyridae.current_position = random_in_box(this->corner1, this->corner2, g);

        //set orientation part
        qt tmp_o = c.ligands[0].rigid.orientation;
        quaternion_increment(tmp_o, random_inside_sphere(g));
        single_lampyridae.current_orientation = tmp_o;

        //init. the array for the number of torsion
        single_lampyridae.current_torsion = new fl[this->torsionSize];
        //single_lampyridae.pbest_torsion=new fl[this->torsionSize];

        for (int x = 0; x < this->torsionSize; x++) //init. all the torsion that the ligand has
        {
            //single_lampyridae.vT[x] = random_fl(-pi, pi, g);
            single_lampyridae.current_torsion[x] = random_fl(-pi, pi, g);
        }

        fireflies.push_back(single_lampyridae);
    }

    firefly::gbest_fit = 1.7976931348623158e+308;
}

void firefly::moveFireflyPosition(int master, int slave, rng &generator)
{
    vec master_pos = fireflies[master].current_position;
    vec slave_pos = fireflies[slave].current_position;
    fl distance_sqr = sqr(master_pos[0] - slave_pos[0]) + sqr(master_pos[1] - slave_pos[1]) + sqr(master_pos[2] - slave_pos[2]);
    fireflies[slave].current_position[0] += beta * (master_pos[0] - slave_pos[0]) * std::exp(gamma * distance_sqr * (-1)) + alpha * random_fl(0, 1, generator);
    fireflies[slave].current_position[1] += beta * (master_pos[1] - slave_pos[1]) * std::exp(gamma * distance_sqr * (-1)) + alpha * random_fl(0, 1, generator);
    fireflies[slave].current_position[2] += beta * (master_pos[2] - slave_pos[2]) * std::exp(gamma * distance_sqr * (-1)) + alpha * random_fl(0, 1, generator);

    if (fireflies[slave].current_position[0] < corner1[0])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);
    if (fireflies[slave].current_position[1] < corner1[1])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);
    if (fireflies[slave].current_position[2] < corner1[2])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);

    if (fireflies[slave].current_position[0] > corner2[0])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);
    if (fireflies[slave].current_position[1] > corner2[1])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);
    if (fireflies[slave].current_position[2] > corner2[2])
        fireflies[slave].current_position = random_in_box(this->corner1, this->corner2, this->g);
}

void firefly::moveFireflyOrientation(int master, int slave, rng &generator)
{
    qt master_ori = fireflies[master].current_orientation;
    qt slave_ori = fireflies[slave].current_orientation;

    //fl distance_sqr = sqr(master_ori.R_component_1() - slave_ori.R_component_1()) + sqr(master_ori.R_component_2() - slave_ori.R_component_2()) + sqr(master_ori.R_component_3() - slave_ori.R_component_3()) + sqr(master_ori.R_component_4() - slave_ori.R_component_4());

    fl distance_sqr = quaternion_difference(master_ori, slave_ori).norm();
    
    quaternion_increment(fireflies[slave].current_orientation, beta * quaternion_to_angle(master_ori - slave_ori) * std::exp(gamma * distance_sqr * (-1)));

    quaternion_increment(fireflies[slave].current_orientation, vec(alpha * random_fl(0, 1, generator), alpha * random_fl(0, 1, generator), alpha * random_fl(0, 1, generator)));
}

void firefly::moveFireflyTorsion(int master, int slave, rng &generator, sz which)
{
    fl master_tor = fireflies[master].current_torsion[which];
    fl slave_tor = fireflies[slave].current_torsion[which];
    fl distance_sqr = sqr(master_tor - slave_tor);
    fireflies[slave].current_torsion[which] += beta * (master_tor - slave_tor) * std::exp(gamma * distance_sqr * (-1)) + alpha * random_fl(0, 1, generator);
}

vec firefly::getCurrentPosition(int i)
{
    return fireflies[i].current_position;
}

qt firefly::getCurrentOrientation(int i)
{
    return fireflies[i].current_orientation;
}

fl firefly::getCurrentTorsion(int i, sz which)
{
    return fireflies[i].current_torsion[which];
}

void firefly::updateGlobalBest(int i)
{
    firefly::gbest_fit = fireflies[i].current_fit;
}

void firefly::updateCurrentFit(int i, double e)
{
    fireflies[i].current_fit = e;
}

double firefly::getCurrentFit(int i)
{
    return fireflies[i].current_fit;
}

void firefly::updateCurrentPosition(int i, vec pos)
{
    fireflies[i].current_position = pos;
}

void firefly::updateCurrentOrientation(int i, qt orientation)
{
    fireflies[i].current_orientation = orientation;
}

void firefly::updateCurrentTorsion(int i, fl torsion, sz which)
{
    fireflies[i].current_torsion[which] = torsion;
}

void firefly::updateGlobalBestFit(double e)
{
    firefly::gbest_fit = e;
}

void firefly::updatePersonalBest(int i, double e)
{
    fireflies[i].pbest_fit = e;
}

double firefly::getPersonalBest(int i)
{
    return fireflies[i].pbest_fit;
}
