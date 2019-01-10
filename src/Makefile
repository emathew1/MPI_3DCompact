
include Makefile.in
CFLAGS += -I$(DECOMP_DIR)

OBJECTS  = Compact8Filter.o MPI_3DCompact.o Utils.o Pade6.o SpongeBC.o CurvilinearCSolver_Core.o CurvilinearCSolver_BC.o CurvilinearInterpolator.o AbstractSingleBlockMesh.o Options.o

default: MPI_3DCompact

MPI_3DCompact:
	@echo
	@echo "///////////////////////////////////////////"
	@echo " Making Pencil Decomposition Components... " 
	@echo "///////////////////////////////////////////"
	@echo
	make -C $(DECOMP_DIR)
	@echo
	@echo "///////////////////////////////////////////"
	@echo " Making Solver Components..." 
	@echo "///////////////////////////////////////////"
	@echo
	make MPI_3DCompact.exe 


MPI_3DCompact.exe: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ -I$(INC) $(OBJECTS) -L$(LIB) -L$(DECOMP_DIR) -lc2decomp $(LIBF) 

MPI_3DCompact.o: MPI_3DCompact.cpp Domain.hpp BC.hpp Utils.hpp Macros.hpp TimeStepping.hpp AbstractCSolver.hpp AbstractRK.hpp CurvilinearCSolver.hpp TVDRK3.hpp AbstractSingleBlockMesh.hpp AlgebraicSingleBlockMesh.hpp Options.hpp 
	$(CC) $(CFLAGS) -c $<

Utils.o: Utils.cpp Domain.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

Pade6.o: Pade6.cpp Pade6.hpp AbstractDerivatives.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp 
	$(CC) $(CFLAGS) -c $<

SpongeBC.o: SpongeBC.cpp SpongeBC.hpp Domain.hpp BC.hpp Macros.hpp IdealGas.hpp AbstractSingleBlockMesh.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

AbstractSingleBlockMesh.o: AbstractSingleBlockMesh.cpp AbstractSingleBlockMesh.hpp AbstractCSolver.hpp Domain.hpp Pade6.hpp Adt.hpp
	$(CC) $(CFLAGS) -c $<

Compact8Filter.o: Compact8Filter.cpp Compact8Filter.hpp AbstractFilter.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp 
	$(CC) $(CFLAGS) -c $<

CurvilinearCSolver_BC.o: CurvilinearCSolver_BC.cpp CurvilinearCSolver.hpp Macros.hpp Utils.hpp SpongeBC.hpp AbstractCSolver.hpp AlgebraicSingleBlockMesh.hpp PngWriter.hpp Adt.hpp 
	$(CC) $(CFLAGS) -c $<

CurvilinearCSolver_Core.o: CurvilinearCSolver_Core.cpp CurvilinearCSolver.hpp Macros.hpp Utils.hpp SpongeBC.hpp AbstractCSolver.hpp AlgebraicSingleBlockMesh.hpp PngWriter.hpp Adt.hpp Pade6.cpp Compact8Filter.hpp
	$(CC) $(CFLAGS) -c $<

CurvilinearInterpolator.o: CurvilinearInterpolator.cpp CurvilinearInterpolator.hpp Macros.hpp Utils.hpp AbstractCSolver.hpp AlgebraicSingleBlockMesh.hpp BC.hpp Domain.hpp
	$(CC) $(CFLAGS) -c $<

Options.o: Options.cpp Options.hpp TimeStepping.hpp BC.hpp
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf   *.o *.exe *.a

superclean:
	make clean -C $(DECOMP_DIR)
	rm -r    *.o *.exe *.a


