 
FORTRAN=gfortran

main:
	# either compile like this
	#${FORTRAN} -c typesDef.f90
	#${FORTRAN} -c C2PRoutines.f90
	#${FORTRAN} -c PDE.f90
	#${FORTRAN} -c Main.f90
	#${FORTRAN} -o con2prim Main.o C2PRoutines.o  typesDef.o PDE.o

	# or like this
	${FORTRAN} -o con2prim typesDef.f90 C2PRoutines.f90 PDE.f90 Main.f90

clean:
	rm -f *.o *.mod con2prim