/*
        PSOVina version 1.0                     Date: 26/11/2014

        Authors: Marcus C. K. Ng  <marcus.ckng@gmail.com>

                 Shirley W. I. Siu <shirleysiu@umac.mo>


        Computational Biology and Bioinformatics Lab

        University of Macau

        http://cbbio.cis.umac.mo

*/
#include "pso.h"
#include "random.h"

    // the vector of degree of freedom
	qt pso::gbest_orientation;
	fl* pso::gbest_torsion;
	vec pso::gbest_position;

	double pso::gbest_fit;
//	double pso::gbest_fit_1;

	pso::pso(int num_birds,double w,double c1,double c2,const vec corner1,const vec corner2, rng& g,conf& c)
	{

		  sz torsionSize = c.ligands[0].torsions.size();
		  this->w = w;
		  this->c1 = c1;
		  this->c2 = c2;
		  this->number = num_birds;
		  this->g = g;
		  this->corner1[0] = corner1[0];			//minmum
		  this->corner1[1] = corner1[1];
		  this->corner1[2] = corner1[2];
		  this->corner2[0] = corner2[0];			//maximum
		  this->corner2[1] = corner2[1];
		  this->corner2[2] = corner2[2];
		  this->torsionSize = (int)torsionSize;
		  this->R1Max_ = 1;
		  this->R1Min_ = 0;
		  this->R2Max_ = 1;
		  this->R2Min_ = 0;
		  pso::gbest_torsion = new fl[torsionSize];
		  init(g,c);
		
	}
	
	void pso::init(rng &g,conf& c)
	{

		int i;
		for(i=0;i<this->number;i++)
		{
			bird single_bird;

			single_bird.pbest_fit = 1.7976931348623158e+308;
			single_bird.tmp_fit = 1.7976931348623158e+308;
			
			//set position part
			single_bird.velocity = random_in_box(this->corner1,this->corner2,g);
			single_bird.current_position = random_in_box(this->corner1,this->corner2,g);
			
			//set orientation part
			single_bird.vO = random_inside_sphere(g);
			qt tmp_o = c.ligands[0].rigid.orientation;
			quaternion_increment(tmp_o,  random_inside_sphere(g));
			single_bird.current_orientation = tmp_o;
			
			//init. the array for the number of torsion
			single_bird.current_torsion=new fl[this->torsionSize];
			single_bird.vT=new fl[this->torsionSize];
			single_bird.pbest_torsion=new fl[this->torsionSize];
			
			for(int x=0;x<this->torsionSize;x++)						//init. all the torsion that the ligand has
			{
				single_bird.vT[x] = random_fl(-pi, pi, g);
				single_bird.current_torsion[x] = random_fl(-pi, pi, g);
			
			}
			
			particle.push_back(single_bird);					
		}
		
		
		pso::gbest_fit = 1.7976931348623158e+308;

	}
	
	
	void pso::updateVelocity(rng& generator,int i)
	{
			double r1 = random_double(this->R1Min_,this->R1Max_,generator);
			double r2 = random_double(this->R2Min_,this->R2Max_,generator);
			
			//calculate by wv(t) + c1*r1*(pbest-x(t))+c2*r2(gbest-x(t))
			particle[i].velocity[0] = particle[i].velocity[0]*w+c1*r1*(particle[i].pbest_pos[0]-particle[i].current_position[0])+c2*r2*(pso::gbest_position[0]-particle[i].current_position[0]);
			particle[i].velocity[1] = particle[i].velocity[1]*w+c1*r1*(particle[i].pbest_pos[1]-particle[i].current_position[1])+c2*r2*(pso::gbest_position[1]-particle[i].current_position[1]);
			particle[i].velocity[2] = particle[i].velocity[2]*w+c1*r1*(particle[i].pbest_pos[2]-particle[i].current_position[2])+c2*r2*(pso::gbest_position[2]-particle[i].current_position[2]);
				
	}
	
	void pso::updateVelocityO(rng& generator,int i)
	{
	
		    qt p1 = particle[i].pbest_orientation-particle[i].current_orientation;
		    qt p2 = pso::gbest_orientation-particle[i].current_orientation;
			double r1 = random_double(this->R1Min_,this->R1Max_,generator);
			double r2 = random_double(this->R2Min_,this->R2Max_,generator);
			//calculate by wv(t) + c1*r1*(pbest-x(t))+c2*r2(gbest-x(t))
			particle[i].vO = particle[i].vO*this->w+c1*r1*quaternion_to_angle(p1)+c2*r2*quaternion_to_angle(p2);
			
	}
	
	void pso::updateVelocityT(rng& generator,int i,sz which)
	{
			double r1 = random_double(this->R1Min_,this->R1Max_,generator);
			double r2 = random_double(this->R2Min_,this->R2Max_,generator);
			//calculate by wv(t) + c1*r1*(pbest-x(t))+c2*r2(gbest-x(t))
			particle[i].vT[which] = particle[i].vT[which]*this->w+c1*r1*(particle[i].pbest_torsion[which]-particle[i].current_torsion[which])+c2*r2*(pso::gbest_torsion[which]-particle[i].current_torsion[which]);
	}
	
	
	void pso::computeNewPositions(int i)
	{

			particle[i].current_position[0] = particle[i].current_position[0] + particle[i].velocity[0];
			particle[i].current_position[1] = particle[i].current_position[1] + particle[i].velocity[1];
			particle[i].current_position[2] = particle[i].current_position[2] + particle[i].velocity[2];
			
			//give a random position, if outside the search box
			if(particle[i].current_position[0] < corner1[0])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
			if(particle[i].current_position[1] < corner1[1])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
			if(particle[i].current_position[2] < corner1[2])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
			
			if(particle[i].current_position[0] > corner2[0])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
			if(particle[i].current_position[1] > corner2[1])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);
			if(particle[i].current_position[2] > corner2[2])
				particle[i].current_position = random_in_box(this->corner1,this->corner2,this->g);

	}

    void pso::computeNewOrientation(int i)
	{
			vec tmp_v = particle[i].vO;
			quaternion_increment(particle[i].current_orientation, tmp_v);
	}
	void pso::computeNewTorsion(int i,rng& generator,sz which)
	{
			particle[i].current_torsion[which] = particle[i].current_torsion[which] + particle[i].vT[which];

			//if(particle[i].current_torsion[which] > pi)
				//particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
			//else if(particle[i].current_torsion[which] < -pi)
				//particle[i].current_torsion[which] = random_fl(-pi, pi, this->g);
	}
	
	
	void pso::updatePersonalBest(int i,double e)
	{
		particle[i].pbest_fit = e;
	}

	
	void pso::updateGlobalBest(int i)
	{
		pso::gbest_fit = particle[i].pbest_fit;
	}
	
	void pso::updateGlobalBest_1(fl en)
	{
		pso::gbest_fit = en;
	}

	
	double pso::getPersonalBest(int i)
	{
		return particle[i].pbest_fit;
	}

	void pso::updateBestPosition(int i,vec pos)
	{
		particle[i].pbest_pos = pos;
		
	}
	
	void pso::updateBestOrientation(int i, qt orientation)
	{
		particle[i].pbest_orientation = orientation;
	}
	
	void pso::updateBestTorsion(int i, fl torsion,sz which)
	{
		particle[i].pbest_torsion[which] = torsion;
	}

	
	vec pso::getCurrentPosition(int i)
	{
		return particle[i].current_position;
	}
	
	qt pso::getCurrentOrientation(int i)
	{
		return particle[i].current_orientation;
	}
	
	fl pso::getCurrentTorsion(int i,sz which)
	{
		return particle[i].current_torsion[which];
	}

