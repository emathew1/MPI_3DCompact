#include "C2Decomp.hpp"

void C2Decomp::writeOne(int ipencil, double *var, string filename){

    MPI_Offset disp, filesize;
    MPI_File fh;
    MPI_Datatype data_type, new_type;

    int sizes[3], subsizes[3], starts[3];
    
    sizes[0] = decompMain.xsz[0];
    sizes[1] = decompMain.ysz[1];
    sizes[2] = decompMain.zsz[2];

    if(ipencil == 0){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.xsz[ip];
	    starts[ip] = decompMain.xst[ip]-1;
	}
    }else if(ipencil == 1){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.ysz[ip];
	    starts[ip] = decompMain.yst[ip]-1;
	}
    }else if(ipencil == 2){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.zsz[ip];
	    starts[ip] = decompMain.zst[ip]-1;
	}
    }

    data_type = MPI_DOUBLE;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);
    MPI_Type_commit(&new_type);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    filesize = 0;
    MPI_File_set_size(fh, filesize);

    disp = 0;
    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);
    
    MPI_File_write_all(fh, var, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    MPI_Type_free(&new_type);

}

void C2Decomp::readOne(int ipencil, double *var, string filename){
    
    MPI_Offset disp;
    MPI_File fh;
    MPI_Datatype data_type, new_type;
    
    int sizes[3], subsizes[3], starts[3];

    sizes[0] = decompMain.xsz[0];
    sizes[1] = decompMain.ysz[1];
    sizes[2] = decompMain.zsz[2];

    if(ipencil == 0){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.xsz[ip];
	    starts[ip] = decompMain.xst[ip]-1;
	}
    }else if(ipencil == 1){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.ysz[ip];
	    starts[ip] = decompMain.yst[ip]-1;
	}
    }else if(ipencil == 2){
	for(int ip = 0; ip < 3; ip++){
	    subsizes[ip] = decompMain.zsz[ip];
	    starts[ip] = decompMain.zst[ip]-1;
	}
    }

    data_type = MPI_DOUBLE;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);
    MPI_Type_commit(&new_type);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

    disp = 0;
    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);

    MPI_File_read_all(fh, var, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);

    MPI_Type_free(&new_type);
}

void C2Decomp::writeVar(MPI_File &fh, MPI_Offset &disp, int ipencil, double *var){

    int sizes[3], subsizes[3], starts[3];
    MPI_Datatype data_type, new_type;

    sizes[0] = decompMain.xsz[0];
    sizes[1] = decompMain.ysz[1];
    sizes[2] = decompMain.zsz[2];

    if(ipencil == 0){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.xsz[ip];
            starts[ip] = decompMain.xst[ip]-1;
        }
    }else if(ipencil == 1){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.ysz[ip];
            starts[ip] = decompMain.yst[ip]-1;
        }
    }else if(ipencil == 2){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.zsz[ip];
            starts[ip] = decompMain.zst[ip]-1;
        }
    }

    data_type = MPI_DOUBLE;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);
    MPI_Type_commit(&new_type);

    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);

    MPI_File_write_all(fh, var, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);
    
    MPI_Type_free(&new_type);

    disp += sizes[0]*sizes[1]*sizes[2]*myTypeBytes;
    
}

void C2Decomp::readVar(MPI_File &fh, MPI_Offset &disp, int ipencil, double *var){

    int sizes[3], subsizes[3], starts[3];
    MPI_Datatype data_type, new_type;

    sizes[0] = decompMain.xsz[0];
    sizes[1] = decompMain.ysz[1];
    sizes[2] = decompMain.zsz[2];

    if(ipencil == 0){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.xsz[ip];
            starts[ip] = decompMain.xst[ip]-1;
        }
    }else if(ipencil == 1){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.ysz[ip];
            starts[ip] = decompMain.yst[ip]-1;
        }
    }else if(ipencil == 2){
        for(int ip = 0; ip < 3; ip++){
            subsizes[ip] = decompMain.zsz[ip];
            starts[ip] = decompMain.zst[ip]-1;
        }
    }

    data_type = MPI_DOUBLE;

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);
    MPI_Type_commit(&new_type);

    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);

    MPI_File_read_all(fh, var, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);

    MPI_Type_free(&new_type);

    disp += sizes[0]*sizes[1]*sizes[2]*myTypeBytes;

}

void C2Decomp::writeScalar(MPI_File &fh, MPI_Offset &disp, int n, double *var){

    int m;

    MPI_Datatype data_type;
    data_type = MPI_DOUBLE;

    MPI_File_set_view(fh, disp, data_type, data_type, "native", MPI_INFO_NULL);
    
    if(nRank == 0){
	m = n;
    }else{
	m = 0;
    }

    MPI_File_write_all(fh, var, m, data_type, MPI_STATUS_IGNORE);
    
    disp += n*myTypeBytes;
}

void C2Decomp::readScalar(MPI_File &fh, MPI_Offset &disp, int n, double *var){

    MPI_Datatype data_type;
    data_type = MPI_DOUBLE;
   
    MPI_File_set_view(fh, disp, data_type, data_type, "native", MPI_INFO_NULL);

    MPI_File_read_all(fh, var, n, data_type, MPI_STATUS_IGNORE);

    disp += n*myTypeBytes;
}

void C2Decomp::writePlane(int ipencil, double *var, int iplane, int n, string filename){

    double *wk, *wk2;
    double *wk2d = NULL;

    int sizes[3], subsizes[3], starts[3];

    MPI_Datatype data_type, new_type;
    MPI_File fh;
    MPI_Offset filesize, disp;

    data_type = MPI_DOUBLE;

    if(iplane == 0){
	bool allocFlag = 0;
	if(ipencil == 0){
	    wk = var;
	}else if(ipencil == 1){
	    allocX(wk); allocFlag = 1;
	    transposeY2X(var, wk);
	}else if(ipencil == 2){
	    allocX(wk); allocFlag = 1;
	    allocY(wk2);
	    transposeZ2Y(var,wk2);
	    transposeY2X(wk2,wk);
	    deallocXYZ(wk2);
	}
	
	wk2d = new double[decompMain.xsz[1]*decompMain.xsz[2]];
	
	for(int kp = 0; kp < decompMain.xsz[2]; kp++){
	    for(int jp = 0; jp < decompMain.xsz[1]; jp++){
		int ii2 = kp*decompMain.xsz[1] + jp;
		int ii3 = kp*decompMain.xsz[1]*decompMain.xsz[0] + jp*decompMain.xsz[0] + n;
		wk2d[ii2] = wk[ii3];
	    }
	}

	sizes[0] = 1;
	sizes[1] = decompMain.ysz[1];
	sizes[2] = decompMain.zsz[2];
	subsizes[0] = 1;
	subsizes[1] = decompMain.xsz[1];
	subsizes[2] = decompMain.xsz[2];
	starts[0] = 0;
	starts[1] = decompMain.xst[1]-1;
	starts[2] = decompMain.xst[2]-1;
	
	if(allocFlag == 1){
	    deallocXYZ(wk);
	}
    }else if(iplane == 1){
	bool allocFlag = 0;
	if(ipencil == 0){
	    allocY(wk); allocFlag = 1;
	    transposeX2Y(var, wk);
	}else if(ipencil == 1){
	    wk = var;
	}else if(ipencil == 2){
	    allocY(wk); allocFlag = 1;
	    transposeZ2Y(var,wk);
	}
	
	wk2d = new double[decompMain.ysz[0]*decompMain.ysz[2]];
	
	for(int kp = 0; kp < decompMain.ysz[2]; kp++){
	    for(int ip = 0; ip < decompMain.ysz[0]; ip++){
		int ii2 = kp*decompMain.ysz[0] + ip;
		int ii3 = kp*decompMain.ysz[1]*decompMain.ysz[0] + n*decompMain.ysz[0] + ip;
		wk2d[ii2] = wk[ii3];
	    }
	}

	sizes[0] = decompMain.xsz[0];
	sizes[1] = 1;
	sizes[2] = decompMain.zsz[2];
	subsizes[0] = decompMain.ysz[0];
	subsizes[1] = 1;
	subsizes[2] = decompMain.ysz[2];
	starts[0] = decompMain.yst[0]-1;
	starts[1] = 0;
	starts[2] = decompMain.yst[2]-1;
	
	if(allocFlag == 1){
	    deallocXYZ(wk);
	}

    }else if(iplane == 2){
	bool allocFlag = 0;
	if(ipencil == 0){
	    allocZ(wk); allocFlag = 1;
	    allocY(wk2); 
	    transposeX2Y(var, wk2);
	    transposeY2Z(wk2, wk);
	    deallocXYZ(wk2);
	}else if(ipencil == 1){
	    allocZ(wk); allocFlag = 1;
	    transposeY2Z(var, wk);
	}else if(ipencil == 2){
	    wk = var;
	}
	
	wk2d = new double[decompMain.zsz[0]*decompMain.zsz[1]];
	
	for(int jp = 0; jp < decompMain.zsz[1]; jp++){
	    for(int ip = 0; ip < decompMain.zsz[0]; ip++){
		int ii2 = jp*decompMain.zsz[0] + ip;
		int ii3 =  n*decompMain.zsz[1]*decompMain.zsz[0] + jp*decompMain.zsz[0] + ip;
		wk2d[ii2] = wk[ii3];
	    }
	}

	sizes[0] = decompMain.xsz[0];
	sizes[1] = decompMain.ysz[1];
	sizes[2] = 1;
	subsizes[0] = decompMain.zsz[0];
	subsizes[1] = decompMain.zsz[1];
	subsizes[2] = 1;
	starts[0] = decompMain.zst[0]-1;
	starts[1] = decompMain.zst[1]-1;
	starts[2] = 0;
	
	if(allocFlag == 1){
	    deallocXYZ(wk);
	}

    }

    MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, data_type, &new_type);

    MPI_Type_commit(&new_type);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    filesize = 0;
    MPI_File_set_size(fh, filesize);
   
    disp = 0;
    MPI_File_set_view(fh, disp, data_type, new_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, wk2d, subsizes[0]*subsizes[1]*subsizes[2], data_type, MPI_STATUS_IGNORE);
 
    MPI_File_close(&fh);
    MPI_Type_free(&new_type);

    delete[] wk2d;


}

void C2Decomp::writeEvery(int ipencil, double *var, int iskip, int jskip, int kskip, string filename, bool from1){

    double *wk, *wk2;
    MPI_Offset filesize, disp;
    MPI_File fh;
    MPI_Datatype data_type, new_type;
    MPI_Comm newcomm;

    int sizes[3], subsizes[3], starts[3];
    int key, color; 
    
    int xsz[3], ysz[3], zsz[3];
    int xst[3], yst[3], zst[3];
    int xen[3], yen[3], zen[3];
    int skip[3];

    data_type = MPI_DOUBLE;
    
    skip[0] = iskip;
    skip[1] = jskip;
    skip[2] = kskip;

}


