
include Makefile.in
CFLAGS += -I$(DECOMP_DIR)

OBJECTS  = Filter.o MPI_3DCompact.o Utils.o Derivatives.o UniformCSolver.o

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
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) -L$(DECOMP_DIR) -lc2decomp $(LIBF) 

MPI_3DCompact.o: MPI_3DCompact.cpp Domain.hpp BC.hpp Utils.hpp Macros.hpp TimeStepping.hpp AbstractCSolver.hpp UniformCSolver.hpp AbstractRK.hpp TVDRK3.hpp 
	$(CC) $(CFLAGS) -c $<

Utils.o: Utils.cpp Domain.hpp Utils.hpp
	$(CC) $(CFLAGS) -c $<

Derivatives.o: Derivatives.cpp Derivatives.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp 
	$(CC) $(CFLAGS) -c $<

Filter.o: Filter.cpp Filter.hpp Domain.hpp BC.hpp Utils.hpp Macros.hpp 
	$(CC) $(CFLAGS) -c $<

UniformCSolver.o: UniformCSolver.cpp Macros.hpp Utils.hpp SpongeBC.hpp AbstractCSolver.hpp
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf   *.o *.exe *.a

superclean:
	make clean -C $(DECOMP_DIR)
	rm -r    *.o *.exe *.a


