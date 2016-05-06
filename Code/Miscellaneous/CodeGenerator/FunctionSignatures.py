#!/bin/env python
##
#
# Collection of all function signatures
# 

import re

# set via Backend
m_precision = ''


def setPrecision(i_precision):
    global m_precision
    m_precision = i_precision


#
# signatures of the space-time predictor
#
def getPicardLoopSignature(i_nDim):
    # choose function signature prototype
    if(i_nDim==2):
        l_functionSignature = "template<void PDEFlux2d(const DATATYPE* const Q, DATATYPE* f, DATATYPE* g)>\n"              \
                              "void kernels::aderdg::optimised::picardLoop( \n"                                            \
                              "  DATATYPE* restrict lqh, \n"                                                               \
                              "  DATATYPE* restrict lFh, \n"                                                               \
                              "  const DATATYPE* restrict const luh, \n"                                                   \
                              "  const tarch::la::Vector<DIMENSIONS,DATATYPE> &dx,\n"                                      \
                              "  const DATATYPE dt \n"                                                                     \
                              ")"
    elif(i_nDim==3):
        l_functionSignature = "template<void PDEFlux3d(const DATATYPE* const Q, DATATYPE* f, DATATYPE* g, DATATYPE* h)>\n" \
                              "void kernels::aderdg::optimised::picardLoop( \n"                                            \
                              "  DATATYPE* restrict lqi, \n"                                                               \
                              "  DATATYPE* restrict lFh, \n"                                                               \
                              "  const DATATYPE* restrict const luh, \n"                                                   \
                              "  const tarch::la::Vector<DIMENSIONS,DATATYPE> &dx,\n"                                      \
                              "  const DATATYPE dt \n"                                                                     \
                              ")"
    else:
        l_functionSignature = ""
        print("FunctionSignatures.getPicardLoopSignature(): nDim not supported")


    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getPicardLoopSignature(): precision not supported")
    
    
    return l_functionSignature


def getPredictorSignature():
    # function signature prototype
    l_functionSignature = "void kernels::aderdg::optimised::predictor( \n"       \
                          "  DATATYPE* restrict lqhi, \n"                        \
                          "  DATATYPE* restrict lFhi, \n"                        \
                          "  const DATATYPE* restrict const lqh \n"              \
                          "  const DATATYPE* restrict const lFh \n"              \
                          ")"

    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getPredictorSignature(): precision not supported")

    return l_functionSignature


def getExtrapolatorSignature():
    # function signature prototype
    l_functionSignature = "void kernels::aderdg::optimised::extrapolator( \n"   \
                          "  DATATYPE* restrict lQbnd, \n"                      \
                          "  DATATYPE* restrict lFbnd, \n"                      \
                          "  const DATATYPE* restrict const lqhi, \n"           \
                          "  const DATATYPE* restrict const lFhi \n"            \
                          ")"

    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getExtrapolatorSignature(): precision not supported")

    return l_functionSignature


#
# signature of the volume integral
#
def getVolumeIntegralSignature():
    # function signature prototype
    l_functionSignature = "void kernels::aderdg::optimised::volumeIntegral( \n" \
                          "  DATATYPE* restrict lduh, \n"                      \
                          "  const DATATYPE* restrict const lFhi, \n"          \
                          "  const tarch::la::Vector<DIMENSIONS,DATATYPE> &dx\n"\
                          ")"
    
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getVolumeIntegralSignature(): precision not supported")                          

    return l_functionSignature



#
# signature of the surface integral
#
def getSurfaceIntegralSignature():
    # function signature prototype
    l_functionSignature = "void kernels::aderdg::optimised::surfaceIntegral( \n" \
                          "  DATATYPE* restrict lduh, \n"                       \
                          "  const DATATYPE* restrict const lFbnd, \n"          \
                          "  const tarch::la::Vector<DIMENSIONS,DATATYPE> &dx\n" \
                          ")"
    
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getSurfaceIntegralSignature(): precision not supported")                          

    return l_functionSignature



#
# signature of the element update
#
def getSolutionUpdateSignature():
    # function signature prototype
    l_functionSignature = "void kernels::aderdg::optimised::solutionUpdate( \n" \
                          "  DATATYPE* restrict luh, \n"                       \
                          "  const DATATYPE* restrict const lduh, \n"          \
                          "  const DATATYPE dt\n"                               \
                          ")" 
    
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getSolutionUpdateSignature(): precision not supported")
                                 
    return l_functionSignature



def getSolutionAdjustmentSignature():
    # function signature prototype
    l_functionSignature = "template <void PDESolutionAdjustment(const DATATYPE* const x,const DATATYPE J_w,const DATATYPE t,const DATATYPE dt,DATATYPE* Q)> \n" \
                          "void solutionAdjustment(\n" \
                          "DATATYPE* luh,\n" \
                          "const tarch::la::Vector<DIMENSIONS,DATATYPE>& center,\n" \
                          "const tarch::la::Vector<DIMENSIONS,DATATYPE>& dx,\n" \
                          "const DATATYPE t,\n" \
                          "const DATATYPE dt\n" \
                          ")"
                          
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getSolutionAdjustmentSignature(): precision not supported")
                                 
    return l_functionSignature


#
# signature of the Riemann solver
#
def getRiemannSolverSignature():
    # function signature prototype
    l_functionSignature = "template <void PDEEigenvalues(const DATATYPE* const Q, const int normalNonZero, DATATYPE* lambda)>\n" \
                          "void kernels::aderdg::optimised::riemannSolver( \n"                                                   \
                          "  DATATYPE* restrict lFbndL,\n"                                                                          \
                          "  DATATYPE* restrict lFbndR,\n"                                                                          \
                          "  const DATATYPE* restrict const lQbndL,\n"                                                              \
                          "  const DATATYPE* restrict const lQbndR,\n"                                                              \
                          "  const DATATYPE dt,\n"                                                                               \
                          "  const int normalNonZero\n"                                                                          \
                          ")" 
                          
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getRiemannSolverSignature(): precision not supported")
                                 
    return l_functionSignature    


#
# signature of the initial field
#
def getInitialConditionSignature():
    # function signature prototype:
    l_functionSignature = "template <void PDEInitialValues(const DATATYPE* const x,DATATYPE* Q)>\n" \
                          "void initialCondition(\n"                                                \
                          "  DATATYPE* restrict luh,\n"                                             \
                          "  const tarch::la::Vector<DIMENSIONS,DATATYPE>& center,\n"               \
                          "  const tarch::la::Vector<DIMENSIONS,DATATYPE>& dx\n"                    \
                          ")"
    
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getInitialConditionSignature(): precision not supported")
                                 
    return l_functionSignature


#
# signature of the time step computation
# 
def getStableTimeStepSizeSignature():
    # function signature prototype:
    l_functionSignature = "template <void PDEEigenvalues(const DATATYPE* const Q, const int normalNonZero, DATATYPE* lambda)>\n" \
                          "DATATYPE stableTimeStepSize(\n"                                                                         \
                          "  const DATATYPE* restrict const luh,\n"                                                              \
                          "  const tarch::la::Vector<DIMENSIONS,DATATYPE>& dx\n"                                                 \
                          ")"
    
    # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                     
    if(m_precision=='SP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'float', l_functionSignature)
    elif(m_precision=='DP'):
        l_functionSignature = re.sub(r'\bDATATYPE\b', 'double', l_functionSignature)
    else:
        print("FunctionSignatures.getStableTimeStepSizeSignature(): precision not supported")
                                 
    return l_functionSignature    


