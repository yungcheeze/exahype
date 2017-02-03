# *********************************************************************************************
# README
# *********************************************************************************************
#
# Available configuration Parameters for ExaHyPE
#
# export variable  |  default-value  |  further values         |  description
#--------------------------------------------------------------------------------------------------
# ARCHITECTURE        CPU               Phi, KNL, HSW             Hardware-platform
# COMPILER            Intel             GNU                       Used compiler (and linker)
# MODE                Release           Debug, Profile, Asserts   Verbosity and Debug level
# SHAREDMEM           None              OMP, TBB                  Shared-memory parallelisation
# DISTRIBUTEDMEM      None              MPI                       Distributed-memory parallelisation
# BOUNDARYCONDITIONS  None              Periodic                  Type of boundary conditions
# *********************************************************************************************
PEANO_KERNEL_PEANO_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./Peano/peano
PEANO_KERNEL_TARCH_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./Peano/tarch
PEANO_TOOLBOX_MULTISCALELINKEDCELL_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./Peano/multiscalelinkedcell
PEANO_TOOLBOX_SHAREDMEMORY_ORACLES_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./Peano/sharedmemoryoracles
PEANO_TOOLBOX_MPI_BLANCING_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./Peano/mpibalancing
EXAHYPE_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./ExaHyPE
PROJECT_PATH=/mnt/beegfs/osorio/ExaHyPE-Engine/./ApplicationExamples/SRMHD_LimitingADERDG
EXECUTABLE=ExaHyPE-MHD 


PROJECT_CFLAGS+=-DALIGNMENT=16
# Several MPI variants face problems with multithreaded MPI. As we run into 
# such issues multiple times, we disable by default multithreaded MPI in ExaHyE. 
# However, feel free to give it a try in your code on your system by disabling 
# this flag. 
PROJECT_CFLAGS+=-DnoMultipleThreadsMayTriggerMPICalls
PROJECT_CFLAGS+=-DDim3
ifeq ($(SHAREDMEM),)
  SHAREDMEM=TBB
endif
ifeq ($(DISTRIBUTEDMEM),)
  DISTRIBUTEDMEM=MPI
endif




-include /mnt/beegfs/osorio/ExaHyPE-Engine/./ExaHyPE/Makefile




all: 
	@echo MHD 
	@echo =================
	@echo An ExaHyPE solver
