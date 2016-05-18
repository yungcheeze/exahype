#!/bin/env python

import Backend



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



def generateScatter(i_architecture: str, i_nVar: int, i_chunkSize: int) -> str:
    l_signature = 'void scatter(double* restrict in_buf, double* restrict out_buf)'
    l_simdWidth = Backend.m_simdWidth['DP'][i_architecture]
    l_iters = int(i_nVar/l_simdWidth)
    l_remainder = i_nVar % l_simdWidth
    l_offset = Backend.getSizeWithPadding(i_nVar)
    # E.g. 9 Vars => 2*simdWidth + 1 => 2*packed + 1*remainder

    l_startAdress = 0
    l_code = ''
    if(i_architecture == 'hsw'):
        # fully packed
        l_code = '__m256d v1, v2, v3, v4, perm1, perm2, perm3, perm4, res1, res2, res3, res4;\n'
        for it in range(0, l_iters):
            # aligned load
            l_code = l_code + 'v1 = _mm256_load_pd(&in_buf['+str(0*l_offset+l_simdWidth*it)+']);\n'
            l_code = l_code + 'v2 = _mm256_load_pd(&in_buf['+str(1*l_offset+l_simdWidth*it)+']);\n'
            l_code = l_code + 'v3 = _mm256_load_pd(&in_buf['+str(2*l_offset+l_simdWidth*it)+']);\n'
            l_code = l_code + 'v4 = _mm256_load_pd(&in_buf['+str(3*l_offset+l_simdWidth*it)+']);\n'
            # vunpckhpd
            l_code = l_code + 'perm1 = _mm256_unpacklo_pd(v1,v2);\n'
            l_code = l_code + 'perm2 = _mm256_unpackhi_pd(v1,v2);\n'
            l_code = l_code + 'perm3 = _mm256_unpacklo_pd(v3,v4);\n'
            l_code = l_code + 'perm4 = _mm256_unpackhi_pd(v3,v4);\n'
            # vperm2f128
            l_code = l_code + 'res1 = _mm256_permute2f128_pd(perm1,perm3,0b00100000);\n'
            l_code = l_code + 'res3 = _mm256_permute2f128_pd(perm1,perm3,0b00110001);\n'
            l_code = l_code + 'res2 = _mm256_permute2f128_pd(perm2,perm4,0b00100000);\n'
            l_code = l_code + 'res4 = _mm256_permute2f128_pd(perm2,perm4,0b00110001);\n'
            # unaligned store. Should later on become _mm256_store_pd()
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res1);\n'
            l_startAdress = l_startAdress + i_chunkSize
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res2);\n'
            l_startAdress = l_startAdress + i_chunkSize
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res3);\n'
            l_startAdress = l_startAdress + i_chunkSize
            l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res4);\n'
            l_startAdress = l_startAdress + i_chunkSize

        # we need some scalar instructions
        if(l_remainder > 0):
            # aligned packed load - in_buf is padded
            l_startRemainder = l_iters * l_simdWidth
            l_code = l_code + 'v1 = _mm256_load_pd(&in_buf['+str(0*l_offset+l_startRemainder)+']);\n'
            l_code = l_code + 'v2 = _mm256_load_pd(&in_buf['+str(1*l_offset+l_startRemainder)+']);\n'
            l_code = l_code + 'v3 = _mm256_load_pd(&in_buf['+str(2*l_offset+l_startRemainder)+']);\n'
            l_code = l_code + 'v4 = _mm256_load_pd(&in_buf['+str(3*l_offset+l_startRemainder)+']);\n'
            # vunpckhpd
            l_code = l_code + 'perm1 = _mm256_unpacklo_pd(v1,v2);\n'
            l_code = l_code + 'perm2 = _mm256_unpackhi_pd(v1,v2);\n'
            l_code = l_code + 'perm3 = _mm256_unpacklo_pd(v3,v4);\n'
            l_code = l_code + 'perm4 = _mm256_unpackhi_pd(v3,v4);\n'
            # vperm2f128
            l_code = l_code + 'res1 = _mm256_permute2f128_pd(perm1,perm3,0b00100000);\n'
            l_code = l_code + 'res3 = _mm256_permute2f128_pd(perm1,perm3,0b00110001);\n'
            l_code = l_code + 'res2 = _mm256_permute2f128_pd(perm2,perm4,0b00100000);\n'
            l_code = l_code + 'res4 = _mm256_permute2f128_pd(perm2,perm4,0b00110001);\n'

            if(l_remainder > 0):
                l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res1);\n'
                l_startAdress = l_startAdress + i_chunkSize
            if(l_remainder > 1):
                l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res2);\n'
                l_startAdress = l_startAdress + i_chunkSize
            if(l_remainder > 2):
                l_code = l_code + '_mm256_storeu_pd(&out_buf['+str(l_startAdress)+'],res3);\n'


    print(l_code)
    return l_signature+' {\n'+l_code + '}\n'