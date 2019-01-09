#ifndef _BCH_
#define _BCH_

#include <iostream>
#include "Macros.hpp"

class BC{

    public:

	enum BCType {PERIODIC_SOLVE, DIRICHLET_SOLVE};
	enum BCKind {PERIODIC,
		     INTERNALLY_PERIODIC, 
		     SPONGE,
		     RECT_CURVILINEARSPONGE,
		     CYL_CURVILINEARSPONGE,
		     SPHERICAL_CURVILINEARSPONGE,
		     ADIABATIC_WALL,
		     CONST_T_WALL,
		     MOVING_ADIABATIC_WALL,
		     INLET};

	BCType bcXType, bcYType, bcZType;
	BCKind bcX0, bcX1, bcY0, bcY1, bcZ0, bcZ1;

	double periodicDisp[3][3];

	BC(BCType bcXType, BCKind bcX0, BCKind bcX1, 
	   BCType bcYType, BCKind bcY0, BCKind bcY1,
	   BCType bcZType, BCKind bcZ0, BCKind bcZ1,
	   bool bcPeriodic[3], int mpiRank, double periodicDisp[3][3]){

	    this->bcXType = bcXType;
	    this->bcYType = bcYType;
	    this->bcZType = bcZType;

	    this->bcX0 = bcX0;
	    this->bcX1 = bcX1;
	    this->bcY0 = bcY0;
	    this->bcY1 = bcY1;
	    this->bcZ0 = bcZ0;
	    this->bcZ1 = bcZ1;

	    FOR_I3{
		FOR_J3{
	            this->periodicDisp[i][j] = periodicDisp[i][j];
	        }
	    }

	    IF_RANK0 std::cout << std::endl;
	    IF_RANK0 std::cout << " > Initializing boundary conditions..." << std::endl;
	
	    IF_RANK0 std::cout << " >     X BOUNDARY CONDITIONS    " << std::endl;
	    if(bcXType == PERIODIC_SOLVE){
		IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[0] = true;
	    }else{
		if(bcX0 == SPONGE || bcX0 == RECT_CURVILINEARSPONGE || bcX0 == CYL_CURVILINEARSPONGE){
		   IF_RANK0 std::cout << " > X0=SPG";
		}else if(bcX0 == ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > X0=ABW";
		}else if(bcX0 == MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > X0=MOW";
		}

		if(bcX1 == SPONGE || bcX1 == RECT_CURVILINEARSPONGE || bcX1 == CYL_CURVILINEARSPONGE){
		   IF_RANK0 std::cout << "----------------SPG=X1" << std::endl;
		}else if(bcX1 == ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------ABW=X1" << std::endl;
		}else if(bcX1 == MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------MOW=X1" << std::endl;
		}

		bcPeriodic[0] = false;
	    }

	    IF_RANK0 std::cout << " >     Y BOUNDARY CONDITIONS    " << std::endl;
	    if(bcYType == PERIODIC_SOLVE){
	        IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[1] = true;
	    }else{
		if(bcY0 == SPONGE || bcY0 == RECT_CURVILINEARSPONGE || bcY0 == CYL_CURVILINEARSPONGE){
		   IF_RANK0 std::cout << " > Y0=SPG";
		}else if(bcY0 == ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > Y0=ABW";
		}else if(bcY0 == MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > Y0=MOW";
		}

		if(bcY1 == SPONGE || bcY1 == RECT_CURVILINEARSPONGE || bcY1 == CYL_CURVILINEARSPONGE){
		   IF_RANK0 std::cout << "----------------SPG=Y1" << std::endl;
		}else if(bcY1 == ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------ABW=Y1" << std::endl;
		}else if(bcY1 == MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------MOW=Y1" << std::endl;
		}

		bcPeriodic[1] = false;
	    }

	    IF_RANK0 std::cout << " >     Z BOUNDARY CONDITIONS    " << std::endl;
	    if(bcZType == PERIODIC_SOLVE){
		IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[2] = true;
	    }else{
		if(bcZ0 == SPONGE || bcZ0 == RECT_CURVILINEARSPONGE || bcZ0 == CYL_CURVILINEARSPONGE){
		    IF_RANK0 std::cout << " > Z0=SPG";
		}else if(bcZ0 == ADIABATIC_WALL){
		    IF_RANK0 std::cout << " > Z0=ABW";
		}else if(bcZ0 == MOVING_ADIABATIC_WALL){
		    IF_RANK0 std::cout << " > Z0=MOW";
		}

		if(bcZ1 == SPONGE || bcZ1 == RECT_CURVILINEARSPONGE || bcZ1 == CYL_CURVILINEARSPONGE){
		    IF_RANK0 std::cout << "----------------SPG=Z1" << std::endl;
		}else if(bcZ1 == ADIABATIC_WALL){
		    IF_RANK0 std::cout << "----------------ABW=Z1" << std::endl;
		}else if(bcZ1 == MOVING_ADIABATIC_WALL){
		    IF_RANK0 std::cout << "----------------MOW=Z1" << std::endl;
		}
		bcPeriodic[2] = false;
	     }

	     IF_RANK0{
		std::cout << std::endl;
	     }	

	
	     //Clear out the displacement values if not periodic
	     if(bcPeriodic[0] == false){
		periodicDisp[0][0] = 0.0;
		periodicDisp[0][1] = 0.0;
		periodicDisp[0][2] = 0.0;
	     }

	     if(bcPeriodic[1] == false){
		periodicDisp[1][0] = 0.0;
		periodicDisp[1][1] = 0.0;
		periodicDisp[1][2] = 0.0;
	     }

	     if(bcPeriodic[2] == false){
		periodicDisp[2][0] = 0.0;
		periodicDisp[2][1] = 0.0;
		periodicDisp[2][2] = 0.0;
	     }

	}

};

#endif
