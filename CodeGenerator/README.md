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
| lFhi | (nDim+1) * (_nDof_, nDof, nDof, **nVar**) | (nDim+1) * (_nDof_, nDof, nDof, **nVarPad**) | lFhi has nDim+1 blocks |
| lQhi | _nDof_, nDof, nDof, **nVar** | _nDof_, nDof, nDof, **nVarPad** | |