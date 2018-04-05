#include "C2Decomp.hpp"

void C2Decomp::updateHalo(double *in, double *&out, int level, int ipencil){

    int xs = 0, ys = 0, zs = 0, xe = 0, ye = 0, ze = 0;
   
    int s1 = 0, s2 = 0, s3 = 0;

    int icount, ilength, ijump;

    MPI_Request requests[4];
    MPI_Status  status[4];

    int tag_e, tag_w, tag_n, tag_s, tag_t, tag_b;

    MPI_Datatype datatype = MPI_DOUBLE;
    MPI_Datatype halo12, halo21, halo31, halo32;
 
    if(ipencil == 0){ //X Pencil
	s1 = decompMain.xsz[0];
	s2 = decompMain.xsz[1];
	s3 = decompMain.xsz[2];

	xs = 1;
	xe = s1;
	ys = 1 - level;
	ye = s2 + level;
	zs = 1 - level;
	ze = s3 + level;
	
    }else if(ipencil == 1){ //Y Pencil
	s1 = decompMain.ysz[0];
	s2 = decompMain.ysz[1];
	s3 = decompMain.ysz[2];

	xs = 1 - level;
	xe = s1 + level;
	ys = 1;
	ye = s2;
	zs = 1 - level;
	ze = s3 + level;

    }else if(ipencil == 2){ //Z Pencil
	s1 = decompMain.zsz[0];
	s2 = decompMain.zsz[1];
	s3 = decompMain.zsz[2];

	xs = 1 - level;
	xe = s1 + level;
	ys = 1 - level;
	ye = s2 + level;
	zs = 1;
	ze = s3;
    }else{
	int errorcode = 10;
	string errorstring = "Invalid ipencil passed to updateHalo";
	decomp2DAbort(errorcode, errorstring); 
    }
   
    int outXsize = xe - xs + 1; 
    int outYsize = ye - ys + 1; 
    int outZsize = ze - zs + 1; 
    int outsize = outXsize*outYsize*outZsize;
    out = new double[outsize];

    for(int kp = 0; kp < s3; kp++){
	for(int jp = 0; jp < s2; jp++){	
	    for(int ip = 0; ip < s1; ip++){
		int in_ii  = kp*s2*s1 + jp*s1 + ip;
		int out_ii = 0;

		//Get the output array index...
		if(ipencil == 0){
		    out_ii = (kp+level)*outYsize*outXsize + (jp+level)*outXsize + ip;
		}else if(ipencil == 1){
		    out_ii = (kp+level)*outYsize*outXsize + jp*outXsize + (ip+level);
		}else if(ipencil == 2){
		    out_ii = kp*outYsize*outXsize + (jp+level)*outXsize + (ip+level);
		}

		out[out_ii] = in[in_ii];
	    }
	}
    }


 
   //If needed, define MPI derived data type to pack halo data, then call MPI send/receive to exchange data
    if(ipencil == 0){

	//East/West data in local memory already

	//North/South 

	tag_s = coord[0];
	if(coord[0] == dims[0]-1 && periodicY){
	    tag_n = 0;
	}else{
	    tag_n = coord[0]+1;
	}

	icount  = s3 + 2*level;
	ilength = level*s1;
	ijump   = s1*(s2+2*level);
	
	MPI_Type_vector(icount, ilength, ijump, datatype, &halo12);
	MPI_Type_commit(&halo12);

	//Receive from south
	double *southRecvLocation, *northRecvLocation;
	double *southSendLocation, *northSendLocation;

	southRecvLocation = &out[(zs+level-1)*outYsize*outXsize + (ys+level-1)*outXsize 	  + xs-1]; 
	northRecvLocation = &out[(zs+level-1)*outYsize*outXsize + (ye+level-level)*outXsize       + xs-1]; 
	southSendLocation = &out[(zs+level-1)*outYsize*outXsize + (ys+level+level-1)*outXsize     + xs-1]; 
	northSendLocation = &out[(zs+level-1)*outYsize*outXsize + (ye+level-level-level)*outXsize + xs-1]; 

	MPI_Irecv(southRecvLocation,  1, halo12, neighbor[0][3], tag_s, DECOMP_2D_COMM_CART_X, &requests[0]);
	MPI_Irecv(northRecvLocation,  1, halo12, neighbor[0][2], tag_n, DECOMP_2D_COMM_CART_X, &requests[1]);
	MPI_Issend(southSendLocation, 1, halo12, neighbor[0][3], tag_s, DECOMP_2D_COMM_CART_X, &requests[2]);
	MPI_Issend(northSendLocation, 1, halo12, neighbor[0][2], tag_n, DECOMP_2D_COMM_CART_X, &requests[3]);
	MPI_Waitall(4, requests, status);
	MPI_Type_free(&halo12);

	//Top/Bottom

	icount = (s1*(s2+2*level))*level;

	tag_b = coord[1];
	if(coord[1] == dims[1]-1 && periodicZ){
	    tag_t = 0;
	}else{
	    tag_t = coord[1] + 1;
	}

	double *bottomRecvLocation, *topRecvLocation;
	double *bottomSendLocation, *topSendLocation;

	bottomRecvLocation = &out[(zs+level-1)*outYsize*outXsize           + (ys+level-1)*outXsize + xs-1]; 
	topRecvLocation    = &out[(ze+level-level)*outYsize*outXsize       + (ys+level-1)*outXsize + xs-1]; 
	bottomSendLocation = &out[(zs+level+level-1)*outYsize*outXsize     + (ys+level-1)*outXsize + xs-1]; 
	topSendLocation    = &out[(ze+level-level-level)*outYsize*outXsize + (ys+level-1)*outXsize + xs-1]; 

        MPI_Irecv(bottomRecvLocation,  icount, datatype, neighbor[0][5], tag_b, DECOMP_2D_COMM_CART_X, &requests[0]);
	MPI_Irecv(topRecvLocation,     icount, datatype, neighbor[0][4], tag_t, DECOMP_2D_COMM_CART_X, &requests[1]);
	MPI_Issend(bottomSendLocation, icount, datatype, neighbor[0][5], tag_b, DECOMP_2D_COMM_CART_X, &requests[2]);
	MPI_Issend(topSendLocation,    icount, datatype, neighbor[0][4], tag_t, DECOMP_2D_COMM_CART_X, &requests[3]);
	MPI_Waitall(4, requests, status);


    }else if(ipencil == 1){

	//East/West
	tag_w = coord[0];
	if(coord[0] == dims[0]-1 && periodicX){
	    tag_e = 0;
	}else{
	    tag_e = coord[0] + 1;
	}

	icount  = s2*(s3+2*level);
	ilength = level;
	ijump   = s1+2*level;

	MPI_Type_vector(icount, ilength, ijump, datatype, &halo21);
	MPI_Type_commit(&halo21);

	double *eastRecvLocation, *westRecvLocation;
	double *eastSendLocation, *westSendLocation;

	westRecvLocation = &out[(zs+level-1)*outYsize*outXsize + (ys-1)*outXsize + xs+level-1]; 
	eastRecvLocation = &out[(zs+level-1)*outYsize*outXsize + (ys-1)*outXsize + xe+level-level]; 
	westSendLocation = &out[(zs+level-1)*outYsize*outXsize + (ys-1)*outXsize + xs+level+level-1]; 
	eastSendLocation = &out[(zs+level-1)*outYsize*outXsize + (ys-1)*outXsize + xe+level-level-level]; 

	MPI_Irecv(westRecvLocation,  1, halo21, neighbor[1][1], tag_w, DECOMP_2D_COMM_CART_Y, &requests[0]);
	MPI_Irecv(eastRecvLocation,  1, halo21, neighbor[1][0], tag_e, DECOMP_2D_COMM_CART_Y, &requests[1]);
	MPI_Issend(westSendLocation, 1, halo21, neighbor[1][1], tag_w, DECOMP_2D_COMM_CART_Y, &requests[2]);
	MPI_Issend(eastSendLocation, 1, halo21, neighbor[1][0], tag_e, DECOMP_2D_COMM_CART_Y, &requests[3]);
	MPI_Waitall(4, requests, status);
	MPI_Type_free(&halo21);


	//North/South
	//all data is already in local memory, no halo exchange

        //Top/Bottom
	//data is continguous in memory, no need to define derived data type

	tag_b = coord[1];
	if(coord[1]==dims[1]-1 && periodicZ){
	    tag_t = 0;
	}else{
	    tag_t = coord[1] + 1;
	}
	icount = (s2*(s1+2*level))*level;

	MPI_Barrier(MPI_COMM_WORLD);

	double *bottomRecvLocation, *topRecvLocation;
	double *bottomSendLocation, *topSendLocation;

	bottomRecvLocation = &out[(zs+level-1)*outYsize*outXsize           + (ys-1)*outXsize + xs+level-1]; 
	topRecvLocation    = &out[(ze+level-level)*outYsize*outXsize       + (ys-1)*outXsize + xs+level-1]; 
	bottomSendLocation = &out[(zs+level+level-1)*outYsize*outXsize     + (ys-1)*outXsize + xs+level-1]; 
	topSendLocation    = &out[(ze+level-level-level)*outYsize*outXsize + (ys-1)*outXsize + xs+level-1]; 

        MPI_Irecv(bottomRecvLocation,  icount, datatype, neighbor[1][5], tag_b, DECOMP_2D_COMM_CART_Y, &requests[0]);
	MPI_Irecv(topRecvLocation,     icount, datatype, neighbor[1][4], tag_t, DECOMP_2D_COMM_CART_Y, &requests[1]);
	MPI_Issend(bottomSendLocation, icount, datatype, neighbor[1][5], tag_b, DECOMP_2D_COMM_CART_Y, &requests[2]);
	MPI_Issend(topSendLocation,    icount, datatype, neighbor[1][4], tag_t, DECOMP_2D_COMM_CART_Y, &requests[3]);
	MPI_Waitall(4, requests, status);

    }else if(ipencil == 2){

	//East/West

        tag_w = coord[0];
	if(coord[0] == dims[0]-1 && periodicX){
	    tag_e = 0;
	}else{
	    tag_e = coord[0]+1;
	}

	icount  = (s2 + 2*level)*s3;
	ilength = level;
	ijump   = s1+2*level;

	MPI_Type_vector(icount, ilength, ijump, datatype, &halo31);
	MPI_Type_commit(&halo31);

	double *eastRecvLocation, *westRecvLocation;
	double *eastSendLocation, *westSendLocation;

	westRecvLocation = &out[(zs-1)*outYsize*outXsize + (ys+level-1)*outXsize + xs+level-1]; 
	eastRecvLocation = &out[(zs-1)*outYsize*outXsize + (ys+level-1)*outXsize + xe+level-level]; 
	westSendLocation = &out[(zs-1)*outYsize*outXsize + (ys+level-1)*outXsize + xs+level+level-1]; 
	eastSendLocation = &out[(zs-1)*outYsize*outXsize + (ys+level-1)*outXsize + xe+level-level-level]; 

	MPI_Irecv(westRecvLocation,  1, halo31, neighbor[2][1], tag_w, DECOMP_2D_COMM_CART_Z, &requests[0]);
	MPI_Irecv(eastRecvLocation,  1, halo31, neighbor[2][0], tag_e, DECOMP_2D_COMM_CART_Z, &requests[1]);
	MPI_Issend(westSendLocation, 1, halo31, neighbor[2][1], tag_w, DECOMP_2D_COMM_CART_Z, &requests[2]);
	MPI_Issend(eastSendLocation, 1, halo31, neighbor[2][0], tag_e, DECOMP_2D_COMM_CART_Z, &requests[3]);
	MPI_Waitall(4, requests, status);
	MPI_Type_free(&halo31);

	//North/South

	tag_s = coord[1];
	if(coord[1] == dims[1]-1 && periodicY){
	    tag_n = 0;
	}else{
	    tag_n = coord[1]+1;
	}

	icount  = s3;
	ilength = level*(s1 + 2*level);
	ijump   = (s1 + 2*level)*(s2 + 2*level);

	MPI_Type_vector(icount, ilength, ijump, datatype, &halo32); 
	MPI_Type_commit(&halo32);

	//Receive from south
	double *southRecvLocation, *northRecvLocation;
	double *southSendLocation, *northSendLocation;

	southRecvLocation = &out[(zs-1)*outYsize*outXsize + (ys+level-1)*outXsize           + xs+level-1]; 
	northRecvLocation = &out[(zs-1)*outYsize*outXsize + (ye+level-level)*outXsize       + xs+level-1]; 
	southSendLocation = &out[(zs-1)*outYsize*outXsize + (ys+level+level-1)*outXsize     + xs+level-1]; 
	northSendLocation = &out[(zs-1)*outYsize*outXsize + (ye+level-level-level)*outXsize + xs+level-1]; 

	MPI_Irecv(southRecvLocation,  1, halo32, neighbor[2][3], tag_s, DECOMP_2D_COMM_CART_Z, &requests[0]);
	MPI_Irecv(northRecvLocation,  1, halo32, neighbor[2][2], tag_n, DECOMP_2D_COMM_CART_Z, &requests[1]);
	MPI_Issend(southSendLocation, 1, halo32, neighbor[2][3], tag_s, DECOMP_2D_COMM_CART_Z, &requests[2]);
	MPI_Issend(northSendLocation, 1, halo32, neighbor[2][2], tag_n, DECOMP_2D_COMM_CART_Z, &requests[3]);
	MPI_Waitall(4, requests, status);
	MPI_Type_free(&halo32);



    }


}
