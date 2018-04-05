#include "C2Decomp.hpp"

//template <class T> const T& max (const T& a, const T& b) {
 // return (a<b)?b:a;     // or: return comp(a,b)?b:a; for version (2)
//}

void C2Decomp::decomp2DInit(int pRow, int pCol){

	int errorcode, ierr, row, col;

	row = 0; col = 0;

	//Get the mpi rank and size
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &nRank);

	if(pRow == 0 && pCol == 0){
	   best2DGrid(nProc, row, col);
	}else{
	   if(nProc != pRow*pCol){
		errorcode = 1;
		string errorstring = "Invalid 2D processor grid - nproc /= p_row*p_col\n";
		decomp2DAbort(errorcode, errorstring);
	   }else{
		row = pRow;
		col = pCol;
	   }
	}

	dims[0] = row;
	dims[1] = col;

	//Set up the cartesian coordinates...
	periodic[0] = (int)periodicY;
	periodic[1] = (int)periodicZ;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &DECOMP_2D_COMM_CART_X);

	periodic[0] = (int)periodicX;
	periodic[1] = (int)periodicZ;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &DECOMP_2D_COMM_CART_Y);

	periodic[0] = (int)periodicX;
	periodic[1] = (int)periodicY;
	ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, 0, &DECOMP_2D_COMM_CART_Z);

	//Get the current rank's coordinate in the Cartesian communicator...
	ierr = MPI_Cart_coords(DECOMP_2D_COMM_CART_X, nRank, 2, coord);

	//Derive communicators defining sub-groups for Alltoall(v)
	int remain[2] = {1, 0};
	ierr = MPI_Cart_sub(DECOMP_2D_COMM_CART_X, remain, &DECOMP_2D_COMM_COL);

	remain[0] = 0; remain[1] = 1;
	ierr = MPI_Cart_sub(DECOMP_2D_COMM_CART_X, remain, &DECOMP_2D_COMM_ROW);

	//////////////////////////
	//This could be a pitfall right here because of the column major to row major conversion...
	//////////////////////////

	
	//gather information for halo-cell support
	initNeighbor();

	decompInfoInit();

	for(int i = 0; i < 3; i++){
	   //minus 1 to get C style zero start indices
	   xStart[i] = decompMain.xst[i]-1;
	   yStart[i] = decompMain.yst[i]-1;
	   zStart[i] = decompMain.zst[i]-1;
	
	   xEnd[i]   = decompMain.xen[i]-1;
	   yEnd[i]   = decompMain.yen[i]-1;
	   zEnd[i]   = decompMain.zen[i]-1;
	
	   xSize[i]  = decompMain.xsz[i];
	   ySize[i]  = decompMain.ysz[i];
	   zSize[i]  = decompMain.zsz[i];
	}

	ierr = MPI_Type_size(realType, &myTypeBytes);

};

void C2Decomp::decomp2DAbort(int errorCode, string msg){
    int ierr;
    if(!nRank){
        std::cout << "2Decopm_c Error - errorCode: " << errorCode << std::endl;
        std::cout << "Error Message: " << msg << std::endl;
    }
    ierr = MPI_Abort(MPI_COMM_WORLD, errorCode);
};

void C2Decomp::initNeighbor(){

	//X-pencil
	neighbor[0][0] = MPI_PROC_NULL;
	neighbor[0][1] = MPI_PROC_NULL;
	MPI_Cart_shift(DECOMP_2D_COMM_CART_X, 0, 1, &neighbor[0][3], &neighbor[0][2]);
	MPI_Cart_shift(DECOMP_2D_COMM_CART_X, 1, 1, &neighbor[0][5], &neighbor[0][4]); 

	//Y-pencil
	MPI_Cart_shift(DECOMP_2D_COMM_CART_Y, 0, 1, &neighbor[1][1], &neighbor[1][0]);
	neighbor[1][2] = MPI_PROC_NULL;
	neighbor[1][3] = MPI_PROC_NULL;
	MPI_Cart_shift(DECOMP_2D_COMM_CART_Y, 1, 1, &neighbor[1][5], &neighbor[1][4]);

	//Z-pencil
	MPI_Cart_shift(DECOMP_2D_COMM_CART_Z, 0, 1, &neighbor[2][1], &neighbor[2][0]);
	MPI_Cart_shift(DECOMP_2D_COMM_CART_Z, 1, 1, &neighbor[2][3], &neighbor[2][2]);
	neighbor[2][4] = MPI_PROC_NULL;
	neighbor[2][5] = MPI_PROC_NULL;

};

void C2Decomp::decompInfoInit(){

    int bufSize, nx, ny, nz, errorcode;
    nx = nxGlobal; ny = nyGlobal; nz = nzGlobal;

    //verify the global size can actually be distributed as pencils
    if(nx < dims[0] || ny < dims[0] || ny < dims[1] || nz < dims[1]){
	errorcode = 6;
	string msg = "Invalid 2D processor grid. \n Make sure that min(nx, ny) > p_row and min(ny, nz) >= p_col.";
	decomp2DAbort(errorcode, msg);
    }

    if(nx%dims[0] == 0 && ny%dims[0] == 0 && ny%dims[1] == 0 && nz%dims[1] == 0){
	decompMain.even = true;
    }else{
	decompMain.even = false;
    }

    decompMain.x1dist = new int[dims[0]];
    decompMain.y1dist = new int[dims[0]];
    decompMain.y2dist = new int[dims[1]];
    decompMain.z2dist = new int[dims[1]];

    getDist();

    int pdim[3] = {0, 1, 2};
    partition(nx, ny, nz, pdim, decompMain.xst, decompMain.xen, decompMain.xsz);

    pdim[0] = 1; pdim[1] = 0; pdim[2] = 2;
    partition(nx, ny, nz, pdim, decompMain.yst, decompMain.yen, decompMain.ysz);
 
    pdim[0] = 1; pdim[1] = 2; pdim[2] = 0;
    partition(nx, ny, nz, pdim, decompMain.zst, decompMain.zen, decompMain.zsz);

    decompMain.x1cnts = new int[dims[0]];
    decompMain.y1cnts = new int[dims[0]];
    decompMain.y2cnts = new int[dims[1]];
    decompMain.z2cnts = new int[dims[1]];

    decompMain.x1disp = new int[dims[0]];  
    decompMain.y1disp = new int[dims[0]];  
    decompMain.y2disp = new int[dims[1]];  
    decompMain.z2disp = new int[dims[1]];  

    prepareBuffer(&decompMain);

    int size1 = decompMain.xsz[0]*decompMain.xsz[1]*decompMain.xsz[2];
    int size2 = decompMain.ysz[0]*decompMain.ysz[1]*decompMain.xsz[2];
    int size3 = decompMain.zsz[0]*decompMain.zsz[1]*decompMain.zsz[2];

    bufSize = max(size2, size3);
    bufSize = max(size1, bufSize);
    
    if(bufSize > decompBufSize){
	if(work1_r != NULL){
	    delete[] work1_r;
	    work1_r = NULL;
	}

	if(work2_r != NULL){
	    delete[] work2_r;
	    work2_r = NULL;
	}

	work1_r = new double[bufSize];
	work2_r = new double[bufSize];
	decompBufSize = bufSize;
    }

};

void C2Decomp::getDist(){

    int nx, ny, nz;
    nx = nxGlobal; ny = nyGlobal; nz = nzGlobal;

    int *st1 = new int[dims[0]];
    int *en1 = new int[dims[0]];

    distribute(nx, dims[0], st1, en1, decompMain.x1dist);
    distribute(ny, dims[0], st1, en1, decompMain.y1dist);

    delete[] st1;
    delete[] en1;

    int *st2 = new int[dims[1]];
    int *en2 = new int[dims[1]];

    distribute(ny, dims[1], st2, en2, decompMain.y2dist); 
    distribute(nz, dims[1], st2, en2, decompMain.z2dist);

    delete[] st2;
    delete[] en2; 

}

void C2Decomp::distribute(int data1, int proc, int *st, int *en, int *sz){

    int size1, nl, nu;

    size1 = data1/proc;
    nu = data1 - size1*proc;
    nl = proc - nu;

    st[0] = 1;
    sz[0] = size1;
    en[0] = size1;

    for(int i = 1; i < nl; i++){
	st[i] = st[i-1] + size1;
	sz[i] = size1;
	en[i] = en[i-1] + size1;
    }    

    size1 = size1 + 1;

    for(int i = nl; i < proc; i++){
	st[i] = en[i-1] + 1;
	sz[i] = size1;
	en[i] = en[i-1] + size1;
    }

    en[proc-1] = data1;
    sz[proc-1] = data1 - st[proc-1]+1;
};

void C2Decomp::partition(int nx, int ny, int nz, int *pdim, int *lstart, int *lend, int *lsize){

    int gsize;

    for(int i = 0; i < 3; i++){
	
	if(i == 0){
	    gsize = nx;
	}else if(i == 1){
	    gsize = ny;
	}else if(i == 2){
	    gsize = nz;
	}

	if(pdim[i] == 0){
	    lstart[i] = 1;
	    lend[i]   = gsize;
	    lsize[i]  = gsize;

	}else if(pdim[i] == 1){
    	    int *st, *en, *sz;
	    st = new int[dims[0]];
	    en = new int[dims[0]];
	    sz = new int[dims[0]];
	    
	    distribute(gsize, dims[0], st, en, sz);

   	    lstart[i] = st[coord[0]];
	    lend[i]   = en[coord[0]];
	    lsize[i]  = sz[coord[0]];

	    delete[] st;
	    delete[] en;
	    delete[] sz;
	}else if(pdim[i] == 2){

    	    int *st, *en, *sz;
	    st = new int[dims[1]];
	    en = new int[dims[1]];
	    sz = new int[dims[1]];
	    
	    distribute(gsize, dims[1], st, en, sz);

   	    lstart[i] = st[coord[1]];
	    lend[i]   = en[coord[1]];
	    lsize[i]  = sz[coord[1]];

	    delete[] st;
	    delete[] en;
	    delete[] sz;
	}

    }
    

};

void C2Decomp::prepareBuffer(DecompInfo *dii){

    //MPI_Alltoallv buffer info

    for(int i = 0; i < dims[0]; i++){
	dii->x1cnts[i] = dii->x1dist[i]*dii->xsz[1]*dii->xsz[2];
	dii->y1cnts[i] = dii->y1dist[i]*dii->ysz[0]*dii->ysz[2];
	if(i == 0){
	    dii->x1disp[i] = 0;
	    dii->y1disp[i] = 0;
	}else{
	    dii->x1disp[i] = dii->x1disp[i-1] + dii->x1cnts[i-1];
	    dii->y1disp[i] = dii->y1disp[i-1] + dii->y1cnts[i-1];
	}
    }

    for(int i = 0; i < dims[1]; i++){
	dii->y2cnts[i] = dii->ysz[0]*dii->y2dist[i]*dii->ysz[2];
	dii->z2cnts[i] = dii->zsz[0]*dii->zsz[1]*dii->z2dist[i];
	if(i==0){
	    dii->y2disp[i] = 0;
	    dii->z2disp[i] = 0;
	}else{
	    dii->y2disp[i] = dii->y2disp[i-1] + dii->y2cnts[i-1];
	    dii->z2disp[i] = dii->z2disp[i-1] + dii->z2cnts[i-1];
	}

    }

    dii->x1count = dii->x1dist[dims[0]-1]*dii->y1dist[dims[0]-1]*dii->xsz[2];
    dii->y1count = dii->x1count;

    dii->y2count = dii->y2dist[dims[1]-1]*dii->z2dist[dims[1]-1]*dii->zsz[0];
    dii->z2count = dii->y2count;


}

void C2Decomp::decompInfoFinalize(){

    decompBufSize = 0; 

    delete[] decompMain.x1dist;
    delete[] decompMain.y1dist;
    delete[] decompMain.y2dist;
    delete[] decompMain.z2dist;

    delete[] decompMain.x1cnts;
    delete[] decompMain.y1cnts;
    delete[] decompMain.y2cnts;
    delete[] decompMain.z2cnts;

    delete[] decompMain.x1disp;
    delete[] decompMain.y1disp;
    delete[] decompMain.y2disp;
    delete[] decompMain.z2disp;

    delete[] work1_r;
    delete[] work2_r;


    decompMain.x1dist = NULL;
    decompMain.y1dist = NULL;
    decompMain.y2dist = NULL;
    decompMain.z2dist = NULL;

    decompMain.x1cnts = NULL;
    decompMain.y1cnts = NULL;
    decompMain.y2cnts = NULL;
    decompMain.z2cnts = NULL;

    decompMain.x1disp = NULL;
    decompMain.y1disp = NULL;
    decompMain.y2disp = NULL;
    decompMain.z2disp = NULL;

    work1_r = NULL;
    work2_r = NULL;

}




