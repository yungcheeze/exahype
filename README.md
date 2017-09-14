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

### CodeGenerator dependencies (optional) ###

The optimised kernels are optional and can be replaced by the generic ones, thus the CodeGenerator's dependencies are also optional.

Python 3 is required to run the CodeGenerator.

Additionally the CodeGenerator requires:

* libxsmm (https://github.com/hfp/libxsmm) to generate the advanced matrix multiplication code using libxsmm's generator
* Jinja2 (https://github.com/pallets/jinja.git) a python3 template engine to generate the optimised kernel
* MarkupSafe (https://github.com/pallets/markupsafe.git), a dependency from Jinja2.

A script is provided to import all the dependencies locally

Quick installation:

    ./ExaHyPE-Engine/CodeGenerator/importDependenciesLocally.sh 

If Jinja2 is provided by your python installation you can safelly edit the value of the import script configuration parameter ``JINJA2_ALREADY_AVAILABLE`` to remove its and MarkupSafe local imports.


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


## Regenerate ExaHyPE's kernel gluecode ##
 
```
java -jar ~/workspace/peano/pdt/pdt.jar --generate-gluecode exahype/exahype.specification exahype ~/workspace/peano/pdt/usrtemplates:../Peano/multiscalelinkedcell
```


## Build a new release ##

I assume that the ExaHyPE release repository is checked out to ~/git/ExaHyPE-Release.

0. Ensure that your ExaHyPE Release folder is up to date 
   - Usually, git pull should do 
   - If in doubt, try
     git fetch origin
     git reset --hard origin/master
     
1. Update the guidebook
    - Change into the directory holding your guidebook and build with `make release`.
      The release target builds the PDF without annotations.
    - Copy the PDF over: `cp guidebook.pdf ~/git/ExaHyPE-Release`

2. Build the toolkit
    - Change into Toolkit
      `./build.sh && cp dist/* ~/git/ExaHyPE-Release`

3. Create the two repository images
   - Change into your exahype engine's repository:

     ```
     tar -czhvf ExaHyPE.tar.gz --exclude=.svn --exclude=*.o Peano ExaHyPE LICENSE.txt 
     tar -czvf ExaHyPE-without-Peano.tar.gz --exclude=.svn --exclude=*.o --exclude Peano/peano --exclude Peano/tarch Peano ExaHyPE LICENSE.txt 
     mv *.tar.gz ~/git/ExaHyPE-Release
     ```

4. Copy over the source files
    - Change into your exahype engine's repository:
    - Create directories (only once):
    - `mkdir ~/git/ExaHyPE-Release/ExaHyPE`
    - `mkdir ~/git/ExaHyPE-Release/Peano`
    - `mkdir ~/git/ExaHyPE-Release/Toolkit`
    - Actual copy command:

      ```
cp -R LICENSE.txt ~/git/ExaHyPE-Release
cp -R ExaHyPE ~/git/ExaHyPE-Release
cp -R Peano/mpibalancing ~/git/ExaHyPE-Release/Peano
cp -R Peano/multiscalelinkedcell ~/git/ExaHyPE-Release/Peano
cp -R Peano/sharedmemoryoracles ~/git/ExaHyPE-Release/Peano
cp -R Toolkit/src ~/git/ExaHyPE-Release/Toolkit
cp -R Toolkit/src/Manifest.txt ~/git/ExaHyPE-Release/Toolkit/src
cp -R Toolkit/build.sh ~/git/ExaHyPE-Release/Toolkit
```

5. Cleanup
    - Change into the release directory
    - `find . -name "*.o" -delete`
    - `find . -name "*.class" -delete`

6. Clarify/ensure that a new snapshot of Peano is uploaded to Peano's page

7. Push 

8. Log into http://github.com
    - Change into the repository view and click on the tab releases
    - Create a new release
    - Add the tars to the release (at least the two ExaHyPE tars plus the toolkit. And then probably the demonstrators, too.


x. Demonstrators
  The demonstrators are not part of the above description. If you want to 
  release new demonstrators, too, you have to change into the respective 
  directory and create the demonstrator snapshot manually.
  
