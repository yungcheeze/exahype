 
FORTRAN=gfortran
binary=con2prim-standalone

main:
	${FORTRAN} -o $(binary) Parameters.f90 C2PRoutines.f90  C2P-SRHD.f90  Main.f90

clean:
	rm -f *.o *.mod $(binary)
