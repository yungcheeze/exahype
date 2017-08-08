Dependencies
============

The CodeGenerator requires the template engine Jinja2 (http://jinja.pocoo.org/)

If not installed with your python:

1) Clone the source from the git repository: 
		git clone https://github.com/pallets/jinja.git 
 OR
   unpack the provided archive (release 2.9.6) and rename the directory to jinja: 
		tar -xzf jinja-2.9.6.tar.gz && mv jinja-2.9.6 jinja
		

2) Modify TemplatingUtils.py to use the local version of Jinja2

Codegenerator
=============

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
