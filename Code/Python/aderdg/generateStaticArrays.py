"""
.. module:: generateLookupTable
  :platform: Unix, Windows, Mac
  :synopsis: Generates header files defining static arrays up to a specific order.
.. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>

:synopsis: Generates header files defining static arrays up to a specific order.
"""
import numpy as np
from os import remove
from numpy import linalg
import sys
import os
import errno

from fileWriter import writeVectorToFile, writeMatrixToFile
from aderdg import *

def sizeWithPadding(SIMD_SIZE, vectorSize):
    return SIMD_SIZE * ((vectorSize + (SIMD_SIZE-1)) / SIMD_SIZE);
        
#-----------------------------------------------------------
# main
#-----------------------------------------------------------
maxOrder = 8;
minOrder = 3;
minDim   = 2;
maxDim   = 3;

# width of the vector unit, the export uses only double precision, hence:
# architecture    SIMD_SIZE
#    wsm            2
#    snb            4
#    hsw            4
#    knl            8
SIMD_SIZE = 4

# for a matrix-matrix multiplication C = A *B
# the row count of the matrix A has to be a multiple
# of the SIMD_SIZE. This is achieved through padding
# with zeros. As the system matrices currently
# are on the right side (matrix B) no padding is required.
# This is likely to change in the future. 

outputDirectory = "generatedCode"

# create directory for output files if not existing
try:
    os.makedirs(outputDirectory)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise

# remove all .cpp files (we are in append mode!) 
for fileName in os.listdir(outputDirectory):
    if fileName.endswith(".cpp"):
        remove(outputDirectory+"/"+fileName)
       

# later on there may be several output files, 
# one for each padding scheme
filename=outputDirectory+"/matrices.cpp"

# file setup
with open(filename,"a+") as out:
    out.write("#include \"EulerFlow/dg/DGMatrices.h\" \n\n")
out.close


# we process one order after another
# (1) write the #ifdef construct to the output file
# (2) compute the quadrature points and weights
# (3) compute the system matrices and write them to the output file

order = minOrder
while (order <= maxOrder):    
    dim = minDim
    
   
    # Gauss-Legendre nodes and weights
    x, w = np.polynomial.legendre.leggauss(order+1)
    # map onto [0,1]
    xGPN = 0.5*(x+1)
    wGPN = 0.5*w
    
    weights = np.outer(xGPN, xGPN)  
    
    
    # prepare file
    if order == minOrder:
        with open(filename,"a") as out:
            out.write("#if EXAHYPE_ORDER=="+str(order)+"\n")
    else:
        with open(filename,"a") as out:
            out.write("#elif EXAHYPE_ORDER=="+str(order)+"\n")
    
    
        
    #----------------------------------------------------------------
    # compute matrices and export
    #----------------------------------------------------------------
    Kxi = assembleStiffnessMatrix(xGPN, wGPN, order)
    #                (data, n,      m,       name,  aligned?, namespace, output file)
    writeMatrixToFile(Kxi, order+1, order+1, "Kxi", True, "exahype::dg::", filename)
    #writeMatrixToFile(np.transpose(Kxi), order+1, order+1, "Kxi_transp", "exahype::dg::", filename) 
    
    MM = assembleMassMatrix(xGPN, wGPN, order)
    #writeMatrixToFile(MM, order+1, order+1, "MM", filename)
    
    F0 = assembleTimeFluxMatrixF0(xGPN, order)
    writeVectorToFile(F0, order+1, "F0", filename)
    
    K1 = assembleK1(Kxi, xGPN, order) 
    iK1 = np.linalg.inv(K1)
    writeMatrixToFile(iK1, order+1, order+1, "iK1", True, "exahype::dg::", filename)
    
    FLCoeff, _ = BaseFunc1D(0.0, xGPN, order)
    FRCoeff, _ = BaseFunc1D(1.0, xGPN, order)
    FCoeff = [FLCoeff, FRCoeff]
    writeVectorToFile(FLCoeff, order+1, "FLCoeff", filename)
    writeVectorToFile(FRCoeff, order+1, "FRCoeff", filename)
    writeMatrixToFile(FCoeff, 2, order+1, "FCoeff", True, "exahype::dg::", filename)
    
    dudx = assembleDiscreteDerivativeOperator(MM, Kxi)
    writeMatrixToFile(dudx, order+1, order+1, "dudx", True, "exahype::dg::", filename)
    
    # dim == 2 only
    subOutputMatrix = assembleSubOutputMatrix(xGPN, order, dim)
    writeMatrixToFile(subOutputMatrix, (order+1)**dim, (order+1)**dim, "subOutputMatrix", True, "exahype::dg::", filename)

    
    print("Done with order "+ str(order))
    order+=1

with open(filename,"a") as out:
    out.write("#endif \n")
