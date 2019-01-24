#ifndef _BCH_
#define _BCH_

#include <iostream>
#include "Macros.hpp"
#include "Options.hpp"

class BC{

    public:

	Options::BCType bcXType, bcYType, bcZType;
	Options::BCKind bcX0, bcX1, bcY0, bcY1, bcZ0, bcZ1;

	double periodicDisp[3][3];

	BC(Options *opt, bool bcPeriodic[3]){

	    int mpiRank = opt->mpiRank;

	    this->bcXType = opt->bcXType;
	    this->bcYType = opt->bcYType;
	    this->bcZType = opt->bcZType;

	    this->bcX0 = opt->bcX0;
	    this->bcX1 = opt->bcX1;
	    this->bcY0 = opt->bcY0;
	    this->bcY1 = opt->bcY1;
	    this->bcZ0 = opt->bcZ0;
	    this->bcZ1 = opt->bcZ1;

	    FOR_I3{
		FOR_J3{
	            this->periodicDisp[i][j] = opt->periodicDisp[i][j];
	        }
	    }

	    IF_RANK0 std::cout << std::endl;
	    IF_RANK0 std::cout << " > Initializing boundary conditions..." << std::endl;
	
	    IF_RANK0 std::cout << " >     X BOUNDARY CONDITIONS    " << std::endl;
	    if(bcXType == Options::PERIODIC_SOLVE){
		IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[0] = true;
	    }else{
		if(bcX0 == Options::SPONGE){
		   IF_RANK0 std::cout << " > X0=SPG";
		}else if(bcX0 == Options::ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > X0=ABW";
		}else if(bcX0 == Options::MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > X0=MOW";
		}

		if(bcX1 == Options::SPONGE){
		   IF_RANK0 std::cout << "----------------SPG=X1" << std::endl;
		}else if(bcX1 == Options::ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------ABW=X1" << std::endl;
		}else if(bcX1 == Options::MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------MOW=X1" << std::endl;
		}

		bcPeriodic[0] = false;
	    }

	    IF_RANK0 std::cout << " >     Y BOUNDARY CONDITIONS    " << std::endl;
	    if(bcYType == Options::PERIODIC_SOLVE){
	        IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[1] = true;
	    }else{
		if(bcY0 == Options::SPONGE){
		   IF_RANK0 std::cout << " > Y0=SPG";
		}else if(bcY0 == Options::ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > Y0=ABW";
		}else if(bcY0 == Options::MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << " > Y0=MOW";
		}

		if(bcY1 == Options::SPONGE){
		   IF_RANK0 std::cout << "----------------SPG=Y1" << std::endl;
		}else if(bcY1 == Options::ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------ABW=Y1" << std::endl;
		}else if(bcY1 == Options::MOVING_ADIABATIC_WALL){
		   IF_RANK0 std::cout << "----------------MOW=Y1" << std::endl;
		}

		bcPeriodic[1] = false;
	    }

	    IF_RANK0 std::cout << " >     Z BOUNDARY CONDITIONS    " << std::endl;
	    if(bcZType == Options::PERIODIC_SOLVE){
		IF_RANK0 std::cout << " > ----------PERIODIC---------- " << std::endl;
		bcPeriodic[2] = true;
	    }else{
		if(bcZ0 == Options::SPONGE){
		    IF_RANK0 std::cout << " > Z0=SPG";
		}else if(bcZ0 == Options::ADIABATIC_WALL){
		    IF_RANK0 std::cout << " > Z0=ABW";
		}else if(bcZ0 == Options::MOVING_ADIABATIC_WALL){
		    IF_RANK0 std::cout << " > Z0=MOW";
		}

		if(bcZ1 == Options::SPONGE){
		    IF_RANK0 std::cout << "----------------SPG=Z1" << std::endl;
		}else if(bcZ1 == Options::ADIABATIC_WALL){
		    IF_RANK0 std::cout << "----------------ABW=Z1" << std::endl;
		}else if(bcZ1 == Options::MOVING_ADIABATIC_WALL){
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
