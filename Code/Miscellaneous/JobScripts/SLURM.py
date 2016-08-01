from subprocess import call

def runMultipleSlURMjobs(dimensions, processes, pdegrees, hmaxs, compilers, modes):
  for dimension in dimensions:
    for process in processes:
      for pdegree in pdegrees:
        for hmax in hmaxs:
          for compiler in compilers:
            for mode in modes:
              runSingleSlURMjob(dimension, process, pdegree, hmax, compiler, mode)
            
def runSingleSlURMjob(dimension, process, pdegree, hmax, compiler, mode):
            
  name = 'ExaHyPE_Euler__dimension_' + dimension + '__process_' + `process` + '__pdegree_' + `pdegree` + '__hmax_' + hmax + '__compiler_' + compiler + '__mode_' + mode
  filename = name + '.slurm'
  
  file = open(filename, 'w')
  
  file.write("#!/bin/bash"                                                                                                         + "\n")
  file.write("#SBATCH -o %j." + name + ".out "                                                                                     + "\n")
  file.write("#SBATCH -D /home/hpc/pr63so/gu89tik2/jobs/logs"                                                                   + "\n")
  file.write("#SBATCH -J " + name                                                                                                  + "\n")
  file.write("#SBATCH --get-user-env "                                                                                             + "\n")
  file.write("#SBATCH --clusters=mpp2"                                                                                             + "\n")
  file.write("# alternatively, use mpp1 "                                                                                          + "\n")
  file.write("#SBATCH --ntasks=" + `process`                                                                                       + "\n")
  file.write("# multiples of 28 for mpp2"                                                                                          + "\n")
  file.write("# multiples of 16 for mpp1"                                                                                          + "\n")
  # file.write("#SBATCH --cpus-per-task=28"                                                                                         + "\n")
  file.write("#SBATCH --cpus-per-task=1"                                                                                           + "\n")
  file.write("#SBATCH --mail-type=all "                                                                                            + "\n")
  file.write("#SBATCH --mail-user=varduhn@tum.de "                                                                                 + "\n")
  file.write("#SBATCH --export=NONE "                                                                                              + "\n")
  file.write("#SBATCH --time=01:00:00 "                                                                                            + "\n")
  file.write("source /etc/profile.d/modules.sh"                                                                                    + "\n")
  
  file.write(""                                                                                                                    + "\n")

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
  # file.write("export"                                                                                                            + "\n")

  file.write("cd ~/storage/ExaHyPE"                                                                                                + "\n")
  file.write("git pull"                                                                                                            + "\n")
  file.write("cd"                                                                                                                  + "\n")
  
  file.write(""                                                                                                                    + "\n")
  file.write("ExaHyPEtempDir=$(dd if=/dev/random bs=16 count=1 2>/dev/null | od -An -tx1 | tr -d ' \t\n')"                         + "\n")
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

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*5**********************\""                                                                                   + "\n")
  file.write("cd Toolkit"                                                                                                          + "\n")
  file.write("./build.sh"                                                                                                          + "\n")
  file.write("cd .."                                                                                                               + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*6**********************\""                                                                                   + "\n")
  file.write("cd ApplicationExamples"                                                                                              + "\n")
  file.write("echo \"\" > myUserSpec.exahype"                                                                                      + "\n")
  file.write("echo \"/**                                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \" 2D Euler Flow                                                      \" >> myUserSpec.exahype"                 + "\n")
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
    file.write("echo \"    width                    = 1.0, 1.0                           \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    offset                   = 0.0, 0.0                           \" >> myUserSpec.exahype"                 + "\n")
  else:
    file.write("echo \"    dimension                = 3                                  \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    width                    = 1.0, 1.0, 1.0                      \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    offset                   = 0.0, 0.0, 0.0                      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    end-time                 = 0.4                                  \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end computational-domain                                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"/*                                                                  \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  shared-memory                                                     \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    identifier               = autotuning                           \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    cores                    = 2                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    properties-file          = sharedmemory.properties       		   \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end shared-memory                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"*/                                                                  \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  if process > 1:
    file.write("echo \"  distributed-memory                                                \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    identifier               = static_load_balancing                \" >> myUserSpec.exahype"                 + "\n")
    file.write("echo \"    configure                = {greedy,FCFS}                        \" >> myUserSpec.exahype"                 + "\n")
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
  file.write("echo \"    order              = " + `pdegree` + "                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    maximum-mesh-size  = " + hmax + "                               \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    time-stepping      = global                                     \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    kernel             = generic::fluxes::nonlinear                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    language           = C                                          \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    plot vtk::Cartesian::ascii                                      \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"      variables = 5                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"      time      = 0.0                                               \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"      repeat    = 0.4                                               \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"      output    = ./solution                                        \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"      select    = {all}                                             \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"    end plot                                                        \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"  end solver                                                        \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"                                                                    \" >> myUserSpec.exahype"                 + "\n")
  file.write("echo \"end exahype-project                                                 \" >> myUserSpec.exahype"                 + "\n")
  file.write("cd .."                                                                                                               + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*7**********************\""                                                                                   + "\n")
  file.write("java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/myUserSpec.exahype"                         + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("cd ApplicationExamples/EulerFlow"                                                                                    + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*8**********************\""                                                                                   + "\n")
  file.write("export COMPILER=" + compiler                                                                                         + "\n")
  file.write("export MODE=" + mode                                                                                                 + "\n")
  file.write("make clean"                                                                                                          + "\n")
  file.write("make -j56"                                                                                                           + "\n")
  if process == 1:
    file.write("./ExaHyPE-Euler ../myUserSpec.exahype > " + name + ".$SLURM_JOBID.out"                                             + "\n")
  else:
    file.write("mpiexec -np " + `process` + " ./ExaHyPE-Euler ../myUserSpec.exahype | tee " + name + ".$SLURM_JOBID.out"           + "\n")
    
  file.write("cat " + name + ".$SLURM_JOBID.out"                                                                                   + "\n")
  file.write(""                                                                                                                    + "\n")
  if process > 1:
    file.write("python ../../../Code/Peano/peano/performanceanalysis/merge-log-files.py exahype.log-file " + `process`               + "\n")
  else:
    file.write("cp exahype.log-file " + name + ".merged-exahype.log-file"                                                          + "\n")
  file.write("cp merged-exahype.log-file " + name + ".merged-exahype.log-file"                                                     + "\n")
  if mode == "Profile":
    file.write("module unload gcc"                                                                                                   + "\n")
    file.write("module load python"                                                                                                  + "\n")
    file.write("python ../../../Code/Peano/peano/performanceanalysis/performanceanalysis.py " + name + ".merged-exahype.log-file " + `process` + " 0 "  + "\n")
    file.write(""                                                                                                                    + "\n")
    file.write("scp " + name + ".merged-exahype.log-file* varduhnv@atsccs60.informatik.tu-muenchen.de:~/www-exahype/performance.analysis/" + "\n")
    file.write("ssh varduhnv@atsccs60.informatik.tu-muenchen.de \"chmod -R a+rx ~/www-exahype/performance.analysis\""                + "\n")
  else:
    file.write("python ../../Miscellaneous/JobScripts/scaling/extractPerformanceIndicators.py " + name + ".merged-exahype.log-file > " + name + ".scaling "+ "\n")
    file.write("cp " + name + ".scaling ~/jobs/scaling/$SLURM_JOBID." + name + ".scaling "           + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("cp " + name + ".merged-exahype.log-file ~/jobs/results/$SLURM_JOBID." + name + ".merged-exahype.log-file "           + "\n")
  file.write(""                                                                                                                    + "\n")
  file.write("cd ../.."                                                                                                            + "\n")

  file.write(""                                                                                                                    + "\n")
  file.write("echo \"*9******\"$SLURM_JOBID\"****************\""                                                                   + "\n")


  file.close()
  
  call(["sbatch", filename])
  #print "sbatch " + filename
  call("mv " + filename + " ~/jobs/scripts/" + filename, shell=True)

