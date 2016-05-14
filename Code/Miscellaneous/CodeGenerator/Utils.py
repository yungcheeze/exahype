#!/bin/env python



def generateImmediateDSCAL(i_scalar: float, i_inVectorName: str, i_outVectorName: str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++)\n'  \
             '    '+i_outVectorName+'[i] = '+ str(i_scalar) +' * '+i_inVectorName+'[i];\n'
    return l_code


def generateDSCAL(i_scalarName: str, i_inVectorName: str, i_outVectorName:str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++) \n' \
             '    '+i_outVectorName+'[i] = '+ i_scalarName+' * '+i_inVectorName+'[i];\n'
    return l_code
