#ifndef FIREFLY_H_
#define FIREFLY_H_

#include "common.h"
#include "conf.h"
#include <vector>

class firefly
{	
public:
	struct lampyridae{
		     
		//vec pbest_pos;
		double current_fit;
		vec current_position;
 
		//qt pbest_orientation; 
		qt current_orientation;
		     
		//fl* pbest_torsion;		     
		fl* current_torsion;
		     
		double tmp_fit;  

	};
	
	int torsionSize;
	double beta, gamma, alpha;			// weight, learning coefficient(1&2)
	rng g;

	int number;							//number of firefly
	vec corner1,corner2;				//corners of search box
	
	static vec gbest_position;			//global best of degree of freedom vector
	static qt gbest_orientation;		//global best of degree of freedom vector
	static fl* gbest_torsion;			//global best of degree of freedom vector
	
	static double gbest_fit;			//global best value

	std::vector<lampyridae> fireflies;
	
	
	  double R1Max_;
	  double R1Min_;
	  double R2Max_;
	  double R2Min_;
	
	
	firefly(int,double,double,double,const vec,const vec,rng&,conf&);
	void init(rng&,conf&);
        
	void moveFireflyPosition(int, int, rng&);
    void moveFireflyOrientation(int, int, rng&);
	void moveFireflyTorsion(int, int, rng&, sz);
    
	void updateGlobalBest(int);
    
	double getCurrentFit(int);

	void updateCurrentFit(int, double);

	void updateCurrentPosition(int, vec);
	void updateCurrentOrientation(int, qt);
	void updateCurrentOrientation(int, sz, fl*);

    /*Get current vector*/
	vec getCurrentPosition(int);
    qt getCurrentOrientation(int);
    fl getCurrentTorsion(int,sz);
    
};


#endif /*FIREFLY_H_*/
