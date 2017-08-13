Prerequisites
=============

Python
------

Python 3.3 or more is required.

Dependencies
------------

The CodeGenerator uses the template engine Jinja2 (http://jinja.pocoo.org/)

If it is not provided by your python installation you can instead use the source 
of jinja2 directly by following these steps

1) Either:

* Clone the source from the git repository: ``git clone https://github.com/pallets/jinja.git``
 
* Unpack the provided archive (release 2.9.6) and rename the directory to jinja: 
		``tar -xzf jinja-2.9.6.tar.gz && mv jinja-2.9.6 jinja``
		

2) Modify TemplatingUtils.py to use the local version of Jinja2

From: 
```
# 1. if jinja2 is available with Python
from jinja2 import Template

# 2. if jinja2 sources are locally available
# import sys
# sys.path.insert(0, 'jinja')
# from jinja.jinja2 import Template
```

To 
```
# 1. if jinja2 is available with Python
# from jinja2 import Template

# 2. if jinja2 sources are locally available
import sys
sys.path.insert(0, 'jinja')
from jinja.jinja2 import Template
```


Paths
-----

Every path is relative to the root of the project (inside the ExaHyPE-Engine directory).

The CodeGenerator assumes the following:

* it is located in ``CodeGenerator/``
* the internal ExaHyPE is at ``ExaHyPE/``

If this is not the case, you may need to edit

* the configuration parameters of ``Toolkit/src/eu/exahype/CodeGeneratorHelper.java``
* the configuration parameters of ``CodeGenerator/Driver.py``

The generated code will be put accordingly to the ``pathToOptKernel`` argument starting from the internal ExaHyPe, by default in a subdirectory of ``ExaHyPE/kernels/aderdg/optimised/``.


Codegenerator
=============

To access the help: ``python3 CodeGenerator/Drivers.py -h``

Usage and arguments
-------------------

positional arguments:
*  pathToOptKernel    desired relative path to the generated code (../ExaHyPE/ as root)
*  solverName         name of the user-solver
*  numberOfVariables  the number of quantities
*  order              the order of the approximation polynomial
*  dimension          number of dimensions you want to simulate
*  numerics           linear or nonlinear
*  architecture       the microarchitecture of the target device
*  pathToLibxsmm      where to find your local copy of code generator back end 'https://github.com/hfp/libxsmm'

optional arguments:
*  -h, --help         show this help message and exit
*  --deepProfiling    enable deep-rpofiling (use only with profiler enable)
*  --useFlux          enable flux
*  --useNCP           enable non conservative product
*  --useSource        enable source terms
*  --noTimeAveraging  disable time averaging in the spacetimepredictor (less memory usage, more computation)


Example: ``python3 CodeGenerator/Driver.py kernels/aderdg/optimised/test Euler::MyEulerSolver 5 6 2 nonlinear hsw Libxsmm --useFlux``


Data format and padding
-----------------------

The Codegenerator may use padding when producing architecture specific code, it may also change the index order

Using the C index order convention with index in italic being absent in dim 2


| Array | Generic | Optimised | Note |
| ----- | ------- | --------- | ---- | 
| luh | _nDof_, nDof, nDof, nVar | _nDof_, nDof, nDof, nVar | unchanged |
| lFhbnd & lQhbnd | 2*nDim, **_nDof_, nDof**, nVar | 2*nDim, nVar, **_nDof_, nDof** | in 3D the two nDof dim are padded as one block + nVar and nDofs swap |
| lQhi | _nDof_, nDof, nDof, **nVar** | _nDof_, nDof, nDof, **nVarPad** | |
| lFhi | (nDim+1) * (_nDof_, nDof, nDof, **nVar**) | (nDim+1) * (_nDof_, nDof, nDof, **nVarPad**) | lFhi has nDim+1 blocks |
| LQh (LQi+LQhi_old+rhs+rhs_old or tempSpaceTimeUnknowns) | _nDof_, nDof, nDof, nDof, **nVar** | _nDof_, nDof, nDof, nDof, **nVarPad** | nonlinear case |
| LFh (LFi or tempSpaceTimeFluxUnknowns[0]) | nDim+1, _nDof_, nDof, nDof, nDof, **nVar** | nDim+1, _nDof_, nDof, nDof, nDof, **nVarPad** | nonlinear case, +1 for source |
| gradQ (tempSpaceTimeFluxUnknowns[1]) | _nDof_, nDof, nDof, nDof, nDim, nVar | _nDof_, nDof, nDof, nDof, nDim, nVar | Non padded/reordered since used in user code |
