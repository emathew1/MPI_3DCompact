#include "Options.hpp"

void Options::parseTSTypeFromString(string vmKey, string inString, TimeStepping::TimeSteppingType &currentType){

    if(vm.count(vmKey)){
        if(strcmp(inString.c_str(), "CONST_CFL")==0){
            currentType = TimeStepping::CONST_CFL;
        }else if(strcmp(inString.c_str(), "CONST_DT")==0){
            currentType = TimeStepping::CONST_DT;
        }else{
            cout << " > UNKNOWN TIMESTEPPING TYPE SPECIFIED: " << inString << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
        }
        cout << " > " << vmKey << " = " << inString << endl;
    }else{
        cout << " > " << vmKey << " = " << inString << " not specified " << endl;
        MPI_Abort(MPI_COMM_WORLD, -10);
    }
}

void Options::parseBCTypeFromString(string vmKey, string inString, BC::BCType &currentType){
    if(vm.count(vmKey)){
        if(strcmp(inString.c_str(), "PERIODIC_SOLVE")==0){
            currentType = BC::PERIODIC_SOLVE;
        }else if(strcmp(inString.c_str(), "DIRICHLET_SOLVE")==0){
            currentType = BC::DIRICHLET_SOLVE;
        }else{
            cout << " > UNKNOWN BOUNDARY TYPE SPECIFIED: " << inString << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
        }
        cout << " > " << vmKey << " = " << inString << endl;
    }else{
        cout << " > " << vmKey << " = " << inString << " not specified " << endl;
        MPI_Abort(MPI_COMM_WORLD, -10);
    }
}

void Options::parseBCKindFromString(string vmKey, string inString, BC::BCKind &currentType){

    if(vm.count(vmKey)){
        if(strcmp(inString.c_str(), "INTERNALLY_PERIODIC")==0){
            currentType = BC::INTERNALLY_PERIODIC;
        }else if(strcmp(inString.c_str(), "PERIODIC")==0){
            currentType = BC::PERIODIC;
        }else if(strcmp(inString.c_str(), "ADIABATIC_WALL")==0){
            currentType = BC::ADIABATIC_WALL;
        }else if(strcmp(inString.c_str(), "SPONGE")==0){
            currentType = BC::SPONGE;
        }else if(strcmp(inString.c_str(), "CONST_T_WALL")==0){
            currentType = BC::CONST_T_WALL;
        }else if(strcmp(inString.c_str(), "MOVING_ADIABATIC_WALL")==0){
            currentType = BC::MOVING_ADIABATIC_WALL;
        }else if(strcmp(inString.c_str(), "INLET")==0){
            currentType = BC::INLET;
        }else{
            cout << " > UNKNOWN BOUNDARY TYPE SPECIFIED: " << inString << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
        }

        cout << " > " << vmKey << " = " << inString << endl;

    }else{
        cout << " > " << vmKey << " = " << inString << " not specified" << endl;
        MPI_Abort(MPI_COMM_WORLD, -10);
    }
}

void Options::bcValidation(){

   //do X Direction first
   if(bcXType == BC::PERIODIC_SOLVE){
	
	//If not periodic or internally periodic
	if(!(bcX0 == BC::PERIODIC || bcX0 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcXType PERIODIC_SOLVE must have bcX0 be PERIODIC or INTERNALLY PERIODIC, currently bcX0 = " << bcX0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If not periodic or internally periodic
	if(!(bcX1 == BC::PERIODIC || bcX1 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcXType PERIODIC_SOLVE must have bcX1 be PERIODIC or INTERNALLY PERIODIC, currently bcX1 = " << bcX1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If periodic and internally periodic at opposite ends
	if((bcX0 == BC::PERIODIC && bcX1 == BC::INTERNALLY_PERIODIC) || (bcX0 == BC::INTERNALLY_PERIODIC && bcX1 == BC::PERIODIC)){
	    cout << " > bcX0 and bcX1 have mismatched periodic conditions: bcX0 = " << bcX0_str << ", bcX1 = " << bcX1_str << endl; 
            MPI_Abort(MPI_COMM_WORLD, -10);
	} 
   }

   if(bcXType == BC::DIRICHLET_SOLVE){
	//If trying to use periodic conditions in dirichlet solve mode
	if(bcX0 == BC::PERIODIC || bcX0 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcXType DIRICHLET_SOLVE cannot have bcX0 be PERIODIC or INTERNALLY PERIODIC, currently bcX0 = " << bcX0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If trying to use periodic conditions in dirichlet solve mode
	if(bcX1 == BC::PERIODIC || bcX1 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcXType DIRICHLET_SOLVE cannot have bcX1 be PERIODIC or INTERNALLY PERIODIC, currently bcX1 = " << bcX1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}
   }




   //do Y Direction first
   if(bcYType == BC::PERIODIC_SOLVE){
	
	//If not periodic or internally periodic
	if(!(bcY0 == BC::PERIODIC || bcY0 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcYType PERIODIC_SOLVE must have bcY0 be PERIODIC or INTERNALLY PERIODIC, currently bcY0 = " << bcY0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If not periodic or internally periodic
	if(!(bcY1 == BC::PERIODIC || bcY1 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcYType PERIODIC_SOLVE must have bcY1 be PERIODIC or INTERNALLY PERIODIC, currently bcY1 = " << bcY1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If periodic and internally periodic at opposite ends
	if((bcY0 == BC::PERIODIC && bcY1 == BC::INTERNALLY_PERIODIC) || (bcY0 == BC::INTERNALLY_PERIODIC && bcY1 == BC::PERIODIC)){
	    cout << " > bcY0 and bcY1 have mismatched periodic conditions: bcY0 = " << bcY0_str << ", bcY1 = " << bcY1_str << endl; 
            MPI_Abort(MPI_COMM_WORLD, -10);
	} 
   }

   if(bcYType == BC::DIRICHLET_SOLVE){
	//If trying to use periodic conditions in dirichlet solve mode
	if(bcY0 == BC::PERIODIC || bcY0 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcYType DIRICHLET_SOLVE cannot have bcY0 be PERIODIC or INTERNALLY PERIODIC, currently bcY0 = " << bcY0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If trying to use periodic conditions in dirichlet solve mode
	if(bcY1 == BC::PERIODIC || bcY1 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcYType DIRICHLET_SOLVE cannot have bcY1 be PERIODIC or INTERNALLY PERIODIC, currently bcY1 = " << bcY1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}
   }


   //do Z Direction first
   if(bcZType == BC::PERIODIC_SOLVE){
	
	//If not periodic or internally periodic
	if(!(bcZ0 == BC::PERIODIC || bcZ0 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcZType PERIODIC_SOLVE must have bcZ0 be PERIODIC or INTERNALLY PERIODIC, currently bcZ0 = " << bcZ0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If not periodic or internally periodic
	if(!(bcZ1 == BC::PERIODIC || bcZ1 == BC::INTERNALLY_PERIODIC)){
	    cout << " > bcZType PERIODIC_SOLVE must have bcZ1 be PERIODIC or INTERNALLY PERIODIC, currently bcZ1 = " << bcZ1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If periodic and internally periodic at opposite ends
	if((bcZ0 == BC::PERIODIC && bcZ1 == BC::INTERNALLY_PERIODIC) || (bcZ0 == BC::INTERNALLY_PERIODIC && bcZ1 == BC::PERIODIC)){
	    cout << " > bcZ0 and bcZ1 have mismatched periodic conditions: bcZ0 = " << bcZ0_str << ", bcZ1 = " << bcZ1_str << endl; 
            MPI_Abort(MPI_COMM_WORLD, -10);
	} 
   }

   if(bcZType == BC::DIRICHLET_SOLVE){
	//If trying to use periodic conditions in dirichlet solve mode
	if(bcZ0 == BC::PERIODIC || bcZ0 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcZType DIRICHLET_SOLVE cannot have bcZ0 be PERIODIC or INTERNALLY PERIODIC, currently bcZ0 = " << bcZ0_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}

	//If trying to use periodic conditions in dirichlet solve mode
	if(bcZ1 == BC::PERIODIC || bcZ1 == BC::INTERNALLY_PERIODIC){
	    cout << " > bcZType DIRICHLET_SOLVE cannot have bcZ1 be PERIODIC or INTERNALLY PERIODIC, currently bcZ1 = " << bcZ1_str << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
	}
   }

}

void Options::parseSpongeFromString(string vmKey, string inString, SpongeBC::SpongeKind &spongeKind){

    if(vm.count(vmKey)){
        if(strcmp(inString.c_str(), "RECTILINEAR")==0){
            spongeKind = SpongeBC::RECTILINEAR;
        }else if(strcmp(inString.c_str(), "CYLINDRICAL")==0){
            spongeKind = SpongeBC::CYLINDRICAL;
        }else{
            cout << " > UNKNOWN SPONGE TYPE SPECIFIED: " << inString << endl;
            MPI_Abort(MPI_COMM_WORLD, -10);
        }

        cout << " > " << vmKey << " = " << inString << endl;

    }else{
        cout << " > " << vmKey << " = " << inString << " not specified" << endl;
        MPI_Abort(MPI_COMM_WORLD, -10);
    }
}


