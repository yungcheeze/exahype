Codegenerator
=============

Data format and padding
-----------------------

The Codegenerator may use padding when producing architecture specific code

Using the C index order convention with index in italic being absent in dim 2


| Array | Generic | Optimised |
| ----- | ------- | --------- |
| luh | _nDof_, nDof, nDof, nVar | _nDof_, nDof, nDof, nVar |
| lFhbnd & lQhbnd | 2*nDim, _nDof_, nDof, nVar | 2*nDim, **(_nDof_, nDof)** "padded as one block", nVar |
| lFhi | (nDim+1) * (_nDof_, nDof, nDof, **nVar**) | (nDim+1) * (_nDof_, nDof, nDof, **nVarPad**) |
| lQhi | _nDof_, nDof, nDof, **nVar** | _nDof_, nDof, nDof, **nVarPad** |