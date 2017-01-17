# Codegenerator #

## Data format and padding ##

The Codegenerator may use padding when producing architecture specific code

Using the C index order convention with index in italic being absent in dim 2

{|
|+ The table's caption
! Array
! Generic
! Optimised
|-
! luh
| nDof**, nDof, nDof, nVar || nDof**, nDof, nDof, nVar |
|-
! lFhbnd & lQhbnd
| 2 * nDim, ''nDof'', nDof, nVar || | 2 * nDim, '''(''nDof'', nDof) "padded as one block"''', nVar
|-
! lFhi 
| (nDim+1)*(''nDof'', nDof, nDof, '''nVar''') || (nDim+1)*(''nDof'', nDof, nDof, '''nVarPad''') |
|-
! lQhi
| ''nDof'', nDof, nDof, '''nVar''' || ''nDof'', nDof, nDof, '''nVarPad''' |
|}
