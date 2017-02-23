#include "mpi.inc"
#include "tbb.inc"
#include "chrono.inc"

int main(int argc, char** argv) {
	chrono_example();
	tbb_example();
	mpi_example();
}
