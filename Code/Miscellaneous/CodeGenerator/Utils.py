#!/bin/env python



def generateImmediateDSCAL(i_scalar: float, i_vectorName: str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++)\n'  \
             '    '+i_vectorName+'[i] = '+ str(i_scalar) +' * '+i_vectorName+'[i];\n'
    return l_code


def generateDSCAL(i_scalarName: str, i_vectorName: str, i_vectorSize: int) -> str:
    l_code = '#pragma simd\n'                                \
             '  for(int i=0;i<'+str(i_vectorSize)+';i++) \n' \
             '    '+i_vectorName+'[i] = '+ i_scalarName+' * '+i_vectorName+'[i];\n'
    return l_code
