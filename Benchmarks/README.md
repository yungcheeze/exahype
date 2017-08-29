Benchmarks
----------

This folder provides standardised tests for benchmarking the correctness,
as well as the the shared memory and distributed memory
performance of ExaHyPE's generic and optimised kernels.

For some PDEs, e.g. Euler equations, we provide benchmarks
for the pure ADER-DG (Euler_ADERDG) and FV (Euler_FV) implementation in addition
the benchmarks for the coupled ADER-DG - FV method (Euler).

Each implementation directory contains up to four 
subfolders named convergence,multicore,single-node,
and plenty-nodes.

These subfolders contain either two ExaHyPE specification files
with suffix -no-output and exahype, -output.exahype, or,
in case of the convergence folder a single
file with suffix .exahype.
Additionally, there is a job script, e.g. hamilton.slurm-script or
supermuc.load-leveler to be found in the subfolders.

The specification files and the job script serve as templates
for a generator script named <subfolder>/generate-files.sh.
The latter will generate a bunch of job scripts and matching specification
files depending on the chosen test (convergence,multicore,
single-node,plenty-nodes).
The script is intended to be run from the
project directory via
```
<subfolder>/generate-files.sh
```

Modify parameters
-----------------

After the generation, you can still modify the
generated files by hand.
Jobs are also not submitted automatically
you need to that by yourself.

If you want to change some settings for all
generated files, it makes sense to
do that in the original files, i.e.
the job script and specification files
used as templates.

Lastly, feel free to modify the generator
file <subfolder>/generate-files.sh.
However, make sure that the job scripts
still make sense in this case.
Especially when job array scripts are used,
mistakes are easily introduced. 

Cleaning up
-----------

If you want to clean up all generated files,
simply run the clean script
```
<subfolder>/clean.sh
```

Building executables
--------------------
Next to the job and specification file generator scrips,
you can find the scripts configure-(output|no-output).sh
and build-all-(orders|patch-sizes)-(output|no-output).sh.

The configure scripts will run the toolkit against the the project using
one of the subfolder's specification files as input file.

The build-all scripts will do the above for
different ADER-DG orders of approximation or FV patch
sizes and will further call make on the resulting 
source files.
Shared memory and serial builds are created.

It is of course essential to build the 
executables before submitting any jobs.

Removing the executables has to be done manually.
