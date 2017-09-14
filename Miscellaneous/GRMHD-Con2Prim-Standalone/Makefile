 
FORTRAN=gfortran
binary=con2prim-standalone
SET=-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp -Wall -Wno-tabs -Wno-unused-variable -g3

main:
	g++ -c abort.cpp
	${FORTRAN} -cpp -DDim3 $(SET) -c abort.o Parameters.f90 C2PRoutines.f90 C2P-GRMHD.f90 LocalID.f90 PDE.f90 Main.f90
	${FORTRAN} -o $(binary) *.o -lstdc++ 
#  Ccomp/*.o

clean:
	rm -f *.o *.mod $(binary)
