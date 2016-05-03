"""
.. module:: generateLookupTable
  :platform: Unix, Windows, Mac
  :synopsis: Generates lookup table initialisation code up to a specific order.
.. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>,Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Generates lookup table initialisation code up to a specific order.
"""
import numpy as np
from os import remove
from numpy import linalg
import sys
import os
import errno

from fileWriter import writeMatrixLookupTableInitToFile, writeVectorLookupTableInitToFile
from aderdg import *

#-----------------------------------------------------------
# main
#-----------------------------------------------------------
maxOrder = 9;
minOrder = 0;
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
filename2=outputDirectory+"/gausPoints.cpp"

# file setup
with open(filename,"a+") as out:
    out.write("#include \"EulerFlow/dg/DGMatrices.h\" \n\n")
out.close
with open(filename2,"a+") as out2:
    out2.write("#include \"EulerFlow/dg/DGMatrices.h\" \n\n")
out2.close


# we process one order after another
# (1) write the #ifdef construct to the output file
# (2) compute the quadrature points and weights
# (3) compute the system matrices and write them to the output file

out = open("generatedCode/lookupTableInit.csnippet", "w")
out2 = open("generatedCode/gausPoints.csnippet", "w")

order = minOrder
while (order <= maxOrder):    
    dim = minDim
    out.write("// N=%d\n" % order)
    # Gauss-Legendre nodes and weights
    x, w = np.polynomial.legendre.leggauss(order+1)
    # map onto [0,1]
    xGPN = 0.5*(x+1)
    wGPN = 0.5*w
    
    weights = np.outer(xGPN, xGPN)
	
    text = ""
    for l in range(0,order+1):
        line = "gaussLegendreWeights[%d][%d] = %.16g;\n" % (order,l,wGPN[l])
        text += line
    out2.write(text)

    text = ""
    for l in range(0,order+1):
        line = "gaussLegendreNodes  [%d][%d] = %.16g;\n" % (order,l,xGPN[l])
        text += line
    out2.write(text)
	
    out2.write("\n")
    
    #----------------------------------------------------------------
    # compute matrices and export
    #----------------------------------------------------------------
###############################################################################
# Stiffness matrix
###############################################################################
    Kxi  = assembleStiffnessMatrix(xGPN, wGPN, order)
#    Kxi2 = assembleStiffnessMatrixShorter(xGPN, wGPN, order)
#    
#    for j in range(0, order+1):
#        for i in range(0, order+1):
#            print (Kxi[i][j] - Kxi2[i][j])

    #                (data, n,      m,       name,  aligned?, namespace, output file)
    #writeMatrixToFile(Kxi, order+1, order+1, "Kxi", True, "exahype::dg::", filename)
    writeMatrixLookupTableInitToFile( out,"",Kxi,order+1,order+1,"Kxi[%d]" % order );
#    writeMatrixLookupTableInitToFile( out,"",np.transpose(Kxi),order+1,order+1,"Kxi_transposed[%d]" % order);

###############################################################################
# Mass matrix
###############################################################################
#    MM = assembleMassMatrix(xGPN, wGPN, order)
    #writeMatrixToFile(MM, order+1, order+1, "MM", filename)


###############################################################################
# Time lifting operator (is the same as FLCoeff)
###############################################################################
    F0 = assembleTimeFluxMatrixF0(xGPN, order)
    writeVectorLookupTableInitToFile(out, "", F0, order+1, "F0[%d]" % order)

###############################################################################
# Left-hand side matrix appearing in predictor computation 
###############################################################################
    K1 = assembleK1(Kxi, xGPN, order) 
    iK1 = np.linalg.inv(K1)
    writeMatrixLookupTableInitToFile( out,"",iK1,order+1,order+1,"iK1[%d]" % order)

###############################################################################
# Lifting operators
# (Do not see a reason to store FLCoeff and FRCoeff if we have FCoeff)
# (Do not see a reason to store FLCoeff, FRCoeff, and FCoeff if
# we have the equidistiantGridProjectors)
###############################################################################
    FLCoeff, _ = BaseFunc1d(0.0, xGPN, order) 
    FRCoeff, _ = BaseFunc1d(1.0, xGPN, order)
    FCoeff = [FLCoeff, FRCoeff]

    writeVectorLookupTableInitToFile(out, "", FLCoeff, order+1, "FLCoeff[%d]" % order)
    writeVectorLookupTableInitToFile(out, "", FRCoeff, order+1, "FRCoeff[%d]" % order)
    writeMatrixLookupTableInitToFile(out, "", FCoeff, 2, order+1, "FCoeff[%d]" % order)

###############################################################################
# Equidistant grid projectors
###############################################################################
##    dudx = assembleDiscreteDerivativeOperator(MM, Kxi)
##    writeMatrixToFile(dudx, order+1, order+1, "dudx", True, "exahype::dg::", filename)
    # dim == 2 only
#    subOutputMatrix  = assembleSubOutputMatrix(xGPN, order, dim)
#    subOutputMatrix2 = assembleEquidistantGridProjector2d(xGPN, order, dim)    
#    for j in range(0, order+1):
#        for i in range(0, order+1):
#            print (subOutputMatrix[i][j] - subOutputMatrix2[i][j])

#    writeMatrixLookupTableInitToFile(out, "", subOutputMatrix, (order+1)**dim, (order+1)**dim, "subOutputMatrix[%d]" % order)
    equidistantGridProjector1d = assembleEquidistantGridProjector2d(xGPN, order, dim)
    writeMatrixLookupTableInitToFile(out, "", equidistantGridProjector1d, (order+1), (order+1), "equidistantGridProjector1d[%d]" % order)

###############################################################################
# Fine grid projectors
###############################################################################
    fineGridProjector1d0 = assembleFineGridProjector1d(xGPN, 0, order, dim)
    fineGridProjector1d1 = assembleFineGridProjector1d(xGPN, 1, order, dim)
    fineGridProjector1d2 = assembleFineGridProjector1d(xGPN, 2, order, dim)
    writeMatrixLookupTableInitToFile(out, "", fineGridProjector1d0, (order+1), (order+1), "fineGridProjector[%d][0]" % order)
    writeMatrixLookupTableInitToFile(out, "", fineGridProjector1d1, (order+1), (order+1), "fineGridProjector[%d][1]" % order)
    writeMatrixLookupTableInitToFile(out, "", fineGridProjector1d2, (order+1), (order+1), "fineGridProjector[%d][2]" % order)

    print("Done with order "+ str(order))
    order+=1
out.close()


