# This is the ExaHyPE project #

## Mini installation guide ##

Copy and paste these commands to start with a working ExaHyPE application and compile the demo application _EulerFlow_:

    git clone git@gitlab.lrz.de:exahype/ExaHyPE-Engine.git
    ./ExaHyPE-Engine/Peano/checkout-update-peano.sh
    ./ExaHyPE-Engine/Toolkit/build.sh

Now you are ready to follow compile and run an ExaHyPE application [according to the guidebook](http://www5.in.tum.de/exahype/guidebook.pdf):

    cd ExaHyPE-Engine/
    java -jar Toolkit/dist/ExaHyPE.jar ApplicationExamples/EulerFlow.exahype
    cd ApplicationExamples/EulerFlow/
    
    export COMPILER=gnu
    export TBB_INC=/usr/include/tbb
    export TBB_LIB=/usr/lib/tbb
    make
    ./ExaHyPE-Euler ../EulerFlow.exahype

Look into `RUN.sh` for an alternative, more elaborate way how to setp the installation and it's dependencies. You also might want to use the `exa` tool which makes it super simple to start and use an installation from the scratch: After downloading, just execute

    exa=Miscellaneous/BuildScripts/exa.sh 
    $exa bootstrap
    $exa compile-run EulerFlow

### Libxsmm (optional) ###

Libxsmm's code generator is required by the Toolkit when using application taylored optimised kernels to generate the advanced matrix multiplication code. 
The optimised kernels are optional and can be replaced by the generic ones, thus libxsmm is also optional.

Libxsmm's sources can be found on github: https://github.com/hfp/libxsmm

Python 3 is required to compile libxsmm and run the code generator.

Quick installation:

    cd ExaHyPE-Engine/
    git clone https://github.com/hfp/libxsmm.git Libxsmm
    cd Libxsmm/
    make generator

When using the optimised kernels, the path to libxsmm has to be specified in the application's specification file.

The ExaHyPE code generator requires the python3 modules jinja2 and numpy.

## General remarks ##

* Run tests before you commit
* Document your code with doxygen
* Disable auto-formatting of your IDE or follow the google code-style => `Code/Miscellaneous/.clang-format`. For Eclipse users, there is also https://github.com/wangzw/CppStyle.
* Do not run autoformatters on the DaStGen definition files (`*.def`). This will screw up the `Packed-type: ..` and `Constant: ..` lines.


## Commit guidelines ##

Please, don't commit the following:
    
* Binary files (`*.o, executables, ... `) excluding those necessary for the documentation 
* Output files (`*.vtk, logs, ... `)

Please write good commit messages that document how you changed ExaHyPE.



## Build a new release ##

I assume that the ExaHyPE release repository is checked out to ~/git/ExaHyPE-Release. 

1) Update the guidebook
- Change into the directory holding your guidebook and build the pdf.
- Copy the PDF over:
  cp guidebook.pdf ~/git/ExaHyPE-Release

2) Build the toolkit
- Change into Toolkit
- ./build.sh
- cp dist/* ~/git/ExaHyPE-Release

3) Create the two repository images
- Change into your exahype engine's repository:
- tar -czhvf ExaHyPE.tar.gz --exclude=.svn --exclude=*.o Peano ExaHyPE LICENSE.txt 
- tar -czvf ExaHyPE-without-Peano.tar.gz --exclude=.svn --exclude=*.o --exclude Peano/peano --exclude Peano/tarch Peano ExaHyPE LICENSE.txt 
- mv *.tar.gz ~/git/ExaHyPE-Release

4) Copy over the source files
- Change into your exahype engine's repository:
- Create directories (only once):
- mkdir ~/git/ExaHyPE-Release/ExaHyPE
- mkdir ~/git/ExaHyPE-Release/Peano
- mkdir ~/git/ExaHyPE-Release/Toolkit
- Actual copy command:
cp -R LICENSE.txt ~/git/ExaHyPE-Release
cp -R ExaHyPE ~/git/ExaHyPE-Release
cp -R Peano/mpibalancing ~/git/ExaHyPE-Release/Peano
cp -R Peano/multiscalelinkedcell ~/git/ExaHyPE-Release/Peano
cp -R Peano/sharedmemoryoracles ~/git/ExaHyPE-Release/Peano
cp -R Toolkit/src ~/git/ExaHyPE-Release/Toolkit
cp -R Toolkit/src/Manifest.txt ~/git/ExaHyPE-Release/Toolkit/src
cp -R Toolkit/build.sh ~/git/ExaHyPE-Release/Toolkit

5) Cleanup
Change into the release directory
find . -name "*.o" -delete
find . -name "*.class" -delete

6) Clarify/ensure that a new snapshot of Peano is uploaded to Peano's page

7) Push 

8) Log into http://github.com
- Change into the repository view and click on the tab releases
- Create a new release
- Add the tars to the release (at least the two ExaHyPE tars plus the toolkit. And then probably the demonstrators, too.


