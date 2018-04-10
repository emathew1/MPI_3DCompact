
include ../Makefile.in

BASE_OBJECTS = C2Decomp.o Alloc.o TransposeX2Y.o TransposeY2Z.o TransposeZ2Y.o TransposeY2X.o MemSplitMerge.o IO.o Best2DGrid.o Halo.o
TEST_OBJECTS  = test.o 
TESTMAJORINDEXING_OBJECTS  = test_majorindexing.o 
HALOTEST_OBJECTS  = halo_test.o 

default: libc2decomp.a

tests: test_main.exe test_majorindexing.exe test_halo.exe

TransposeX2Y.o: TransposeX2Y.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeY2Z.o: TransposeY2Z.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeZ2Y.o: TransposeZ2Y.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

TransposeY2X.o: TransposeY2X.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Alloc.o: Alloc.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

MemSplitMerge.o: MemSplitMerge.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

IO.o: IO.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Best2DGrid.o: Best2DGrid.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

Halo.o: Halo.cpp C2Decomp.hpp  
	$(CC) $(CFLAGS) -c $< 

C2Decomp.o: C2Decomp.cpp C2Decomp.hpp
	$(CC) $(CFLAGS) -c $< 

test.o: test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

halo_test.o: halo_test.cpp C2Decomp.hpp 
	$(CC) $(CFLAGS) -c $< 

test_main.exe:  $(TEST_OBJECTS) $(BASE_OBJECTS)
	$(CC) $(CFLAGS) $(TEST_OBJECTS) $(BASE_OBJECTS) -o $@ $(LIBF) 

test_majorindexing.exe:  $(TESTMAJORINDEXING_OBJECTS) $(BASE_OBJECTS)
	$(CC) $(CFLAGS) $(TESTMAJORINDEXING_OBJECTS) $(BASE_OBJECTS) -o $@ $(LIBF) 

test_halo.exe:  $(HALOTEST_OBJECTS) $(BASE_OBJECTS)
	$(CC) $(CFLAGS) $(HALOTEST_OBJECTS) $(BASE_OBJECTS) -o $@ $(LIBF) 



libc2decomp.a: $(BASE_OBJECTS)
	$(AR) $@ $^
	$(RANLIB) $@  



clean: 
	rm -rf   *.o *.exe *.a


