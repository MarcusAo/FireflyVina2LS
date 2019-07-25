#ifndef FIREFLY_H_
#define FIREFLY_H_

#include "common.h"
#include "conf.h"
#include <vector>

class firefly
{
  public:
    struct lampyridae
    {

        //vec pbest_pos;
        double current_fit;
        vec current_position;

        //qt pbest_orientation;
        qt current_orientation;

        //fl* pbest_torsion;
        fl *current_torsion;

        double pbest_fit;
    };

    int torsionSize;
    double beta, gamma, alpha; // weight, learning coefficient(1&2)
    double mu1, mu2; //chaotic

    //Levy flight
    double lbeta;
    double sigma_v;
    double sigma_u;

    rng g;

    int number;           //number of firefly
    vec corner1, corner2; //corners of search box

    static vec gbest_position;   //global best of degree of freedom vector
    static qt gbest_orientation; //global best of degree of freedom vector
    static fl *gbest_torsion;    //global best of degree of freedom vector

    static double gbest_fit; //global best values
    static int gbest_firefly;
    std::vector<lampyridae> fireflies;

    double R1Max_;
    double R1Min_;
    double R2Max_;
    double R2Min_;

    firefly(int, double, double, double, double, double, const vec, const vec, rng &, conf &);
    void init(rng &, conf &);

    void moveFireflyPosition(int, int, rng &);
    void moveFireflyOrientation(int, int, rng &);
    void moveFireflyTorsion(int, int, rng &, sz);

    void moveFireflyPosition1(int, int);
    void moveFireflyOrientation1(int, int);
    void moveFireflyTorsion1(int, int, sz);

    void moveFireflyPositionRandomly(int, rng &);
    void moveFireflyOrientationRandomly(int, rng &);
    void moveFireflyTorsionRandomly(int, rng &, sz);

    void updateGlobalBest(int);
    void updateGlobalBestFit(double);
    void updateGlobalBestFirefly(int);

    double getCurrentFit(int);

    void updateCurrentFit(int, double);

    void updateCurrentPosition(int, vec);
    void updateCurrentOrientation(int, qt);
    void updateCurrentTorsion(int, fl, sz);

    /*Get current vector*/
    vec getCurrentPosition(int);
    qt getCurrentOrientation(int);
    fl getCurrentTorsion(int, sz);

    double getPersonalBest(int);
    void updatePersonalBest(int, double);

    double levy(rng &);
    int sign(double);
};

#endif /*FIREFLY_H_*/
