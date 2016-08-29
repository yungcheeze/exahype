from subprocess import call
import os,binascii,datetime
import time

def runMultipleSlURMjobs(dimensions, processes, threads, h_p_ts, compilers, modes, suffix):

  #a tmp directory, where all logfiles go to
  numberOfSimulations = len(dimensions)*len(processes)*len(threads)*len(h_p_ts)*len(compilers)*len(modes)
  directory = time.strftime('%Y_%m_%d__%H_%M_%S') + "_N" + str(numberOfSimulations) + "_" + suffix
  
  for dimension in dimensions:
    for process in processes:
      for thread in threads:
        for h_p_t in h_p_ts:
          for compiler in compilers:
            for mode in modes:
              runSingleSlURMjob(dimension, process, thread, h_p_t, compiler, mode, directory)


def runSingleSlURMjob(dimension, process, thread, h_p_t, compiler, mode, directory):
            
  name = 'ExaHyPE_Euler__' + dimension + '__n_' + `process` + '__t_' + `thread` + '__p_' + `h_p_t[1]` + '__h_' + h_p_t[0] + '__' + compiler + '__' + mode
  filename = name + '.slurm'
  
  file = open(filename, 'w')
  
  writeJobscriptHeader(file, name, thread, process)
  writeJobscriptBody(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)
  
  file.close()
  
  call(["sbatch", filename])
  #print "sbatch " + filename
  #call(["cat", filename])
  call("mv " + filename + " ~/jobs/scripts/" + filename, shell=True)

def writeJobscriptBody(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):  
  writeJobscript_SetupEnvironment(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscript_PrepareSources(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscript_BuildToolkit(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscript_GenerateUserSpecFile(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)
  
  writeJobscript_RunToolkit(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscript_CompileSources(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscript_FilterOutput(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)
  
  writeJobscriptBody_RunExecutable(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)

  writeJobscriptBody_DoPostProcessing(file, name, dimension, process, thread, h_p_t, compiler, mode, directory)
  
              
def writeJobscriptHeader(file, name, thread, process):
  file.write("#!/bin/bash"                                                                                                         + "\n")
  file.write("#SBATCH -o %j." + name + ".out "                                                                                     + "\n")
  file.write("#SBATCH -D /home/hpc/pr63so/gu89tik2/scratch/logs"                                                                   + "\n")
  file.write("#SBATCH -J " + name                                                                                                  + "\n")
  file.write("#SBATCH --get-user-env "                                                                                             + "\n")
  file.write("#SBATCH --clusters=mpp2"                                                                                             + "\n")
  file.write("# alternatively, use mpp1 "                                                                                          + "\n")
  file.write("#SBATCH --ntasks=" + `process`                                                                                       + "\n")
  file.write("# multiples of 28 for mpp2"                                                                                          + "\n")
  file.write("# multiples of 16 for mpp1"                                                                                          + "\n")
  if thread > 1:
    file.write("#SBATCH --cpus-per-task=28"                                                                                           + "\n")
  else:
    file.write("#SBATCH --cpus-per-task=1"                                                                                           + "\n")
  file.write("#SBATCH --mail-type=all "                                                                                            + "\n")
  file.write("#SBATCH --mail-user=varduhn@tum.de "                                                                                 + "\n")
  file.write("#SBATCH --export=NONE "                                                                                              + "\n")
  file.write("#SBATCH --time=01:00:00 "                                                                                            + "\n")
  file.write("source /etc/profile.d/modules.sh"                                                                                    + "\n")
  
def writeJobscript_SetupEnvironment(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write("echo \"*1**********************\""                                                                                   + "\n")
  file.write("[[ $- == *i* ]] && echo 'Interactive' || echo 'Not interactive'"                                                     + "\n")
  file.write("module load gcc"                                                                                                     + "\n")
  file.write("module load git"                                                                                                     + "\n")
  file.write("module load java"                                                                                                     + "\n")
  file.write("module load tbb"                                                                                                     + "\n")
  file.write("module load subversion"                                                                                                     + "\n")
  if compiler == "GNU":
    file.write("module load mpi.ompi/1.10/gcc"                                                                                     + "\n")
    # file.write("module load mpi.intel/5.1_gcc"                                                                                    + "\n")
  else:
    file.write("module load intel"                                                                                                 + "\n")
    file.write("module load mpi.intel"                                                                                             + "\n")
  
  # file.write(""                                                                                                                  + "\n")
  # file.write("export"                                                                                       + "\n")                     
  
def writeJobscript_PrepareSources(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write("cd ~/storage/ExaHyPE"                                                                                                + "\n")
  file.write("git pull"                                                                                                            + "\n")
  file.write("cd"                                                                                                                  + "\n")
  
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*2**********************\""                                                                                   + "\n")
  file.write("cd $SCRATCH/jobs"                                                                                                    + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*3**********************\""                                                                                   + "\n")
  file.write("mkdir $SLURM_JOBID." + name                                                                                          + "\n")
  file.write("cd $SLURM_JOBID." + name                                                                                             + "\n")
  file.write("cp -r ~/storage/ExaHyPE ."                                                                                           + "\n")
  file.write("cd ExaHyPE/Code"                                                                                                     + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*4**********************\""                                                                                   + "\n")
  file.write("cd Peano"                                                                                                            + "\n")
  file.write("tar xfvz peano.tar.gz"                                                                                               + "\n")
  file.write("git checkout .gitignore"                                                                                             + "\n")
  file.write("cd .."                                                                                                               + "\n")


def writeJobscript_BuildToolkit(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*5**********************\""                                                                                   + "\n")
  file.write("cd Toolkit"                                                                                                          + "\n")
  file.write("./build.sh"                                                                                                          + "\n")
  file.write("cd .."                                                                                                               + "\n")

def writeJobscript_GenerateUserSpecFile(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*6**********************\""                                                                                   + "\n")
  file.write("cd ApplicationExamples"                                                                                              + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"\" > myUserSpec.exahype"                                                                                      + "\n")
  file.write("echo \"/**                                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \" 2D/3D Euler Flow                                                   \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \" A simple project                                                   \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \" */                                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"exahype-project  Euler                                              \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  peano-kernel-path          = ./Peano                              \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  exahype-path               = ./ExaHyPE                            \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  output-directory           = ./ApplicationExamples/EulerFlow      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  architecture               = noarch                               \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  computational-domain                                              \" >> myUserSpec.exahype"                 + "\n")
  if dimension == "2D":
    file.write("echo \"    dimension                = 2                                  \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    width                    = 15.0, 15.0                           \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    offset                   = 0.0, 0.0                           \" >> myUserSpec.exahype"                 + "\n")
  else:
    file.write("echo \"    dimension                = 3                                  \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    width                    = 1.0, 1.0, 1.0                      \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    offset                   = 0.0, 0.0, 0.0                      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    end-time                 = " + str(h_p_t[2]) + "              \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end computational-domain                                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  if thread > 1:
    file.write("echo \"  shared-memory                                                     \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    identifier               = dummy                           \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    cores                    = " + `thread` + "                                    \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    properties-file          = sharedmemory.properties       		   \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"  end shared-memory                                                 \" >> myUserSpec.exahype"                 + "\n")
    
  if process > 1:
    file.write("echo \"  distributed-memory                                                \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    identifier               = static_load_balancing                \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    configure                = {hotspot,fair,ranks_per_node:28}                       \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    buffer-size              = 64                                   \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    timeout                  = 120                                  \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"  end distributed-memory                                            \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  optimisation                                                      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    fuse-algorithmic-steps        = on                              \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    fuse-algorithmic-steps-factor = 0.99                            \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end optimisation                                                  \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  solver ADER-DG MyEulerSolver                                      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    variables          = 5                                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    parameters         = 0                                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    order              = " + `h_p_t[1]` + "                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    maximum-mesh-size  = " + h_p_t[0] + "                               \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    time-stepping      = global                                     \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    kernel             = generic::fluxes::nonlinear                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    language           = C                                          \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"    plot vtk::Cartesian::vertices::ascii                            \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"      variables = 5                                                 \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"      time      = 0.0                                               \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"      repeat    = 0.4                                               \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"      output    = ./solution                                        \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"      select    = {all}                                             \" >> myUserSpec.exahype"                 + "\n")
  # file.write("echo \"    end plot                                                        \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end solver                                                        \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"end exahype-project                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("cd .."                                                                                                               + "\n")

def writeJobscript_RunToolkit(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*7**********************\""                                                                                   + "\n")
  file.write("java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/myUserSpec.exahype"                         + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("cd ApplicationExamples/EulerFlow"                                                                                    + "\n")


def writeJobscript_CompileSources(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*8**********************\""                                                                                   + "\n")
  file.write("export COMPILER=" + compiler                                                                                         + "\n")
  file.write("export MODE=" + mode                                                                                                 + "\n")
  file.write("export EXAHYPE_INITIALDATA=DiffusingGauss"                                                                           + "\n")
  file.write("make clean"                                                                                                          + "\n")
  file.write("make -j56"                                                                                                           + "\n")

def writeJobscript_FilterOutput(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write("echo \"# Level Trace Rank Black or white list entry\"                           >> exahype.log-filter" + "\n")
  file.write("echo \"# (info or debug) (-1 means all ranks)\"                                 >> exahype.log-filter" + "\n")
  file.write("echo \"debug tarch -1 black\"                                                   >> exahype.log-filter" + "\n")
  file.write("echo \"debug peano -1 black\"                                                   >> exahype.log-filter" + "\n")
  file.write("echo \"info tarch -1 black\"                                                    >> exahype.log-filter" + "\n")
  file.write("echo \"info peano -1 black\"                                                    >> exahype.log-filter" + "\n")
  file.write("echo \"info peano::utils::UserInterface -1 white\"                              >> exahype.log-filter" + "\n")
  file.write("echo \"info exahype -1 white\"                                                  >> exahype.log-filter" + "\n")
  if mode != "Profile":
    file.write("echo \"\"                                                                       >> exahype.log-filter" + "\n")
    file.write("echo \"info peano::parallel::SendReceiveBufferAbstractImplementation -1 black\" >> exahype.log-filter" + "\n")
  
def writeJobscriptBody_RunExecutable(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  if process == 1:
    file.write("./ExaHyPE-Euler ../myUserSpec.exahype | tee " + name + ".$SLURM_JOBID.out"                                         + "\n")
  else:
    file.write("mpiexec -np " + `process` + " ./ExaHyPE-Euler ../myUserSpec.exahype | tee " + name + ".$SLURM_JOBID.out"           + "\n")
    

def writeJobscriptBody_DoPostProcessing(file, name, dimension, process, thread, h_p_t, compiler, mode, directory):
  file.write(""                                                                                                                    + "\n")
  if process > 1:
    file.write("python ../../../Code/Peano/peano/performanceanalysis/merge-log-files.py exahype.log-file " + `process`               + "\n")
    file.write("cp merged-exahype.log-file " + name + ".merged-exahype.log-file"                                                     + "\n")
  else:
    file.write("cp exahype.log-file " + name + ".merged-exahype.log-file"                                                          + "\n")
  
  if mode == "Profile":
    file.write("module unload gcc"                                                                                                   + "\n")
    file.write("module load python"                                                                                                  + "\n")
    file.write("python ../../../Code/Peano/peano/performanceanalysis/performanceanalysis.py " + name + ".merged-exahype.log-file " + `process` + " 0 "  + "\n")
    file.write(""                                                                                                                    + "\n")
    file.write("scp " + name + ".merged-exahype.log-file* varduhnv@atsccs60.informatik.tu-muenchen.de:~/www-exahype/performance.analysis/" + "\n")
    file.write("ssh varduhnv@atsccs60.informatik.tu-muenchen.de \"chmod -R a+rx ~/www-exahype/performance.analysis\""                + "\n")
  else:
    file.write("python ../../Miscellaneous/JobScripts/scaling/extractPerformanceIndicators.py " + name + ".merged-exahype.log-file > " + name + ".scaling "+ "\n")
    file.write("mkdir ~/jobs/scaling/" + directory                                                                                 + "\n")
    file.write("cp " + name + ".scaling ~/jobs/scaling/" + directory + "/$SLURM_JOBID." + name + ".scaling "                       + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("mkdir ~/jobs/results/" + directory                                                                                 + "\n")
  file.write("cp " + name + ".merged-exahype.log-file ~/jobs/results/" + directory + "/$SLURM_JOBID." + name + ".merged-exahype.log-file " + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("cd ../.."                                                                                                            + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*9******\"$SLURM_JOBID\"****************\""                                                                   + "\n")  
  
  
