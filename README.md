# FireflyVina2LS                                                 
Marcus M. C. Ao & Shirley W. I. Siu                                   

Computational Biology and Bioinformatics Lab (CBBio)                        
Department of Computer and Information Science (CIS)
University of Macau (UM) 
                                                                       
http://cbbio.cis.umac.mo

RELEASE NOTES: 

- FireflyVina2LS was developed based on the framework of PSOVina2LS. 
  For more information about Vina, please visit http://vina.scripps.edu. 

=======================================================================

0. REQUIRED SOFTWARE
   For successful compilation, please install Boost (version 1.59.0) 
   from http://www.boost.org. For preparing molecules for docking, 
   please install AutoDockTools (ADT) from http://mgltools.scripps.edu.

1. INSTALLATION

   The installation basically follows the installation of AutoDock Vina. 
   The steps are simple:

     a. unpack the files
     b. cd firevina2ls/build/<your-platform>/release
     c. modify Makefile to suit your system setting
     d. type "make" to compile
    
    The binary firevina2ls will be generated at the current directory. You can 
    copy this binary to a directory in your PATH e.g. /usr/local/bin, or add
    the path of the current directory to your PATH.

2. RUNNING FireflyVina2LS

   You can run fireflyvina2ls as the way you run vina but additional four 
   parameters (optional) are used to specify how the firefly algorithm performs 
   searching:

   % <path-to-fireflyvina2ls>/fireflyvina2ls

   Firefly parameters (optional):
       --num_fireflies arg (=16)     Number of fireflies
       --beta arg (=1)               Attractiveness
       --gamma arg (=1)              Absorption coefficient 
       --alpha arg (=0.9)           Randomization parameter 
 
   2LS parameters (optional):
       --Cr arg(=15)                 Roughing condition
       --R arg(=0.1)                 Roughing factor

   For example, docking Kifunensine in the Mannosidase enzyme (PDBID 1ps3 from
   the PDBbind v2012 dataset) using FireflyVina2LS with default firefly parameters in a 
   8-core computer and return the lowest energy prediction:

   % <path-to-AutoDockTools>/prepare_ligand4.py -l 1ps3_ligand.mol2 \
     -o 1ps3_ligand.pdbqt -A 'hydrogens' -U 'nphs_lps_waters'

   % <path-to-AutoDockTools>/prepare_receptor4.py -r 1ps3_protein.pdb \
     -o 1ps3_protein.pdbqt -A 'hydrogens' -U 'nphs_lps_waters' 

   % <path-to-fireflyvina2ls>/fireflyvina2ls \
     --receptor 1ps3_protein.pdbqt --ligand 1ps3_ligand.pdbqt \
     --center_x  31.951 --center_y 65.5053 --center_z 7.63888 \
     --size_x    33.452 --size_y   27.612  --size_z   35.136  \ 
     --num_modes 1      --cpu 8  
   
   If you have more CPUs, you can increase the number of particles, 
   e.g. with availability of 16 CPU, try these options:

   --cpu 16 --num_particles 16

3. DEVELOP FireflyVina2LS

   If you are interested in the source code of FireflyVina2LS for any academic
   purposes, please note that the following files were newly developed   
   in our work or modified based on FireflyVina:
   	src/lib/quasi_newton.cpp
	  src/lib/firefly_mutate.cpp

4. CITATION
   Please cite our paper if you have used FireflyVina2LS.

   H. K. Tai, H. Lin and S. W. I. Siu, "Improving the efficiency of PSOVina for 
   protein-ligand docking by two-stage local search," 2016 IEEE Congress on 
   Evolutionary Computation (CEC), Vancouver, BC, 2016, pp. 770-777.

5. CONTACT US
    
   Developer: Marcus M. C. Ao <mingchi_@hotmail.com>
   Project P.I.: Shirley W. I. Siu <shirleysiu@umac.mo>

   Computational Biology and Bioinformatics Lab, University of Macau
   http://cbbio.cis.umac.mo
   http://www.cis.umac.mo/~shirleysiu

