#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Generates the code for the surface integral
# for a specific configuration
#--------------------------------------------------------------
#
#
import Backend
import FunctionSignatures
import Utils


class SurfaceIntegralGenerator:
    m_config = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filename = "surfaceIntegral.cpp"

    # structure of lFbnd, lQbnd in soa-format
    m_chunkSize       = -1
    m_vectorLength    = -1
    m_startAddr_face1 = -1
    m_startAddr_face2 = -1
    m_startAddr_face3 = -1
    m_startAddr_face4 = -1
    m_startAddr_face5 = -1
    m_startAddr_face6 = -1


    def __init__(self, i_config, i_numerics):
        self.m_config = i_config
        self.m_type   = i_numerics

        # access pattern to face data (soa format)
        # lFbnd = [ lFbnd_1 | lFbnd_2 | lFbnd_3 | lFbnd_4 | lFbnd_5 | lFbnd_6 ]
        # lFbnd_1(nDofx,nDofy,nVar)
        self.m_chunkSize       = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        self.m_vectorLength    = self.m_config['nVar'] * self.m_chunkSize
        self.m_startAddr_face1 = 0 * self.m_vectorLength
        self.m_startAddr_face2 = 1 * self.m_vectorLength
        self.m_startAddr_face3 = 2 * self.m_vectorLength
        self.m_startAddr_face4 = 3 * self.m_vectorLength
        self.m_startAddr_face5 = 4 * self.m_vectorLength
        self.m_startAddr_face6 = 5 * self.m_vectorLength


    def __writeHeader(self):
        l_description = '// Solve the surface integral \n\n'

        l_includeStatement = '#include "string.h"\n'                             \
                             '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/aderdg/optimised/DGMatrices.h"\n'\
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n\n'

        l_functionSignature = FunctionSignatures.getSurfaceIntegralSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()




    def generateCode(self):
        self.__writeHeader()

        if(self.m_type == 'linear'):
            self.__generateLinear()
        else:
            self.__generateNonlinear()

        # write missing closing bracket
        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write('}')
        l_sourceFile.close()


    def __generateLinear(self):
        pass


    def __generateNonlinear(self):
        l_sourceFile = open(self.m_filename, 'a')

        # 2D/3D compatibility
        if(self.m_config['nDim'] == 3):
            kMax = self.m_config['nDof']
        else:
            kMax = 1

        # gcc and icc specify distinct ways to inform the compiler about guaranteed alignment
        # gcc: double* arr_ = (double*) __builtin_assume_aligned(a, ALIGNMENT);
        # icc: __assume_aligned(a, ALIGNMENT);
        # the default gcc on the cluster exhibits a well-known bug in alignment assumptions
        # => we skip gcc here
        # do not query __GNUC__ - icc also defines this
        l_sourceFile.write('#ifdef __INTEL_COMPILER\n'\
                           '  __assume_aligned(kernels::s_m, ALIGNMENT)\n'\
                           '  __assume_aligned(kernels::FRCoeff, ALIGNMENT)\n'\
                           '  __assume_aligned(kernels::FLCoeff, ALIGNMENT)\n'
                           '#endif\n')

        # temporary memory for scaled versions
        paddedDof = Backend.getSizeWithPadding(self.m_config['nDof'])
        l_sourceFile.write('  double FRCoeff_s['+str(paddedDof)+'] __attribute__((aligned(ALIGNMENT)));\n'\
                           '  double FLCoeff_s['+str(paddedDof)+'] __attribute__((aligned(ALIGNMENT)));\n')
        # note that FRCoeff and FRCoeff_s may differ w.r.t. length due to padding

        #-----------------------------
        # x-direction
        #-----------------------------
        l_sourceFile.write('  // x-direction\n')
        l_sourceFile.write('  for(int iVar=0;iVar<'+str(self.m_config['nVar'])+';iVar++) {\n')
        l_sourceFile.write('    for(int jk=0;jk<'+str(self.m_config['nDof']**(self.m_config['nDim']-1))+';jk++) {\n')
        # FRCoeff_s = lFbnd(j,k,iVar,2)*FRCoeff
        # copy unpadded length
        l_dscal = Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face2)+'+iVar*'+str(self.m_chunkSize)+'+jk]',\
                                      'kernels::FRCoeff',\
                                      'FRCoeff_s',\
                                      self.m_config['nDof'])
        l_sourceFile.write(l_dscal)
        # FLCoeff_s = lFbnd(j,k,iVar,1)*FLCoeff
        # copy unpadded length
        l_dscal = Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face1)+'+iVar*'+str(self.m_chunkSize)+'+jk]',\
                                      'kernels::FLCoeff',\
                                      'FLCoeff_s',\
                                      self.m_config['nDof'])
        l_sourceFile.write(l_dscal)
        # (weight(j)*weight(k))/dx(1) * (FRCoeff_s - FLCoeff_s)
        # note loop length with padding
        l_sourceFile.write('#pragma simd\n'\
                           '    for(int it=0;it<'+str(paddedDof)+';it++)\n'\
                           '      kernels::s_m[it] = kernels::weights2[jk]/dx[0]*(FRCoeff_s[it]-FLCoeff_s[it]);\n')
        # lduh(iVar,:,j,k) -= s_m(:)
        # scatter...is it worth to use my scatter operation?
        l_stride = self.m_config['nVar']
        for i in range(0, self.m_config['nDof']):
            l_sourceFile.write('  lduh[jk*'+str(self.m_config['nVar']*self.m_config['nDof'])+\
                                       '+'+str(i*l_stride)+\
                                       '+iVar] -= kernels::s_m['+str(i)+'];\n')
            #l_sourceFile.write('  lduh[jk*20+iVar+0*nVar] = kernels::s_m[0]')
        # close for loops
        l_sourceFile.write('    }\n'\
                           '  }\n')


        #-----------------------------
        # y-direction
        #-----------------------------
        l_sourceFile.write('  // y-direction\n')
        if(self.m_config['nDim']==3):
            l_sourceFile.write('  for(int iVar=0;iVar<'+str(self.m_config['nVar'])+';iVar++) {\n')
            l_sourceFile.write('    for(int k=0;k<'+str(self.m_config['nDof'])+';k++) {\n'\
                               '      for(int i=0;i<'+str(self.m_config['nDof'])+';i++) {\n')
            # FRCoeff_s = lFbnd(j,k,iVar,4)*FRCoeff
            l_sourceFile.write(
                Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face4)+'+iVar*'+str(self.m_chunkSize)+'+k*'+str(self.m_config['nDof'])+'+i]',\
                                    'kernels::FRCoeff',\
                                    'FRCoeff_s',\
                                    self.m_config['nDof']))
            # FLCoeff_s = lFbnd(i,k,iVar,3)*FLCoeff
            l_sourceFile.write(
                Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face3)+'+iVar*'+str(self.m_chunkSize)+'+k*'+str(self.m_config['nDof'])+'+i]',\
                                    'kernels::FLCoeff',\
                                    'FLCoeff_s',\
                                    self.m_config['nDof']))

            # (weight(i)*weight(k))/dx(2) * (FRCoeff_s - FLCoeff_s)
            l_sourceFile.write('        double s = kernels::weights1[k]*kernels::weights1[i]/dx[1];\n')
            l_sourceFile.write('#pragma simd\n'\
                               '        for(int it=0;it<'+str(paddedDof)+';it++)\n'\
                               '          kernels::s_m[it] = s*(FRCoeff_s[it]-FLCoeff_s[it]);\n')
            # lduh(iVar,i,:,k) -= s_m(:)
            l_stride = self.m_config['nVar']*self.m_config['nDof']
            for j in range(0, self.m_config['nDof']):
                l_sourceFile.write('  lduh[k*'+str(self.m_config['nVar']*self.m_config['nDof']**2)+\
                                         '+'+str(j*l_stride)+\
                                         '+i*'+str(self.m_config['nVar'])+\
                                         '+iVar] '\
                                   '-= kernels::s_m['+str(j)+'];\n')

            # close for loops
            l_sourceFile.write('      }\n'\
                               '    }\n'\
                               '  }\n')

        elif(self.m_config['nDim']==2):
            print("SurfaceIntegralGenerator.py: 2D y-direction not yet implemented")

        #-----------------------------
        # z-direction
        #-----------------------------
        if(self.m_config['nDim']==3):
            l_sourceFile.write('  // z-direction\n')
            l_sourceFile.write('  for(int iVar=0;iVar<'+str(self.m_config['nVar'])+';iVar++) {\n')
            l_sourceFile.write('    for(int j=0;j<'+str(self.m_config['nDof'])+';j++) {\n'\
                               '      for(int i=0;i<'+str(self.m_config['nDof'])+';i++) {\n')
            # FRCoeff_s = lFbnd(i,j,iVar,6)*FRCoeff
            l_sourceFile.write(
                Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face6)+'+iVar*'+str(self.m_chunkSize)+'+j*'+str(self.m_config['nDof'])+'+i]',\
                                    'kernels::FRCoeff',\
                                    'FRCoeff_s',\
                                    self.m_config['nDof']))
            # FLCoeff_s = lFbnd(i,j,iVar,5)*FLCoeff
            l_sourceFile.write(
                Utils.generateDSCAL('lFbnd['+str(self.m_startAddr_face5)+'+iVar*'+str(self.m_chunkSize)+'+j*'+str(self.m_config['nDof'])+'+i]',\
                                    'kernels::FLCoeff',\
                                    'FLCoeff_s',\
                                    self.m_config['nDof']))

            # (weight(i)*weight(j))/dx(3) * (FRCoeff_s - FLCoeff_s)
            l_sourceFile.write('        double s = kernels::weights1[i]*kernels::weights1[j]/dx[2];\n')
            l_sourceFile.write('#pragma simd\n'\
                               '        for(int it=0;it<'+str(paddedDof)+';it++)\n'\
                               '          kernels::s_m[it] = s*(FRCoeff_s[it]-FLCoeff_s[it]);\n')
            # lduh(iVar,i,j,:) -= s_m(:)
            l_stride = self.m_config['nVar']*(self.m_config['nDof']**2)
            for k in range(0, self.m_config['nDof']):
                l_sourceFile.write('  lduh['+str(k*l_stride)+\
                                         '+j*'+str(self.m_config['nVar']*self.m_config['nDof'])+\
                                         '+i*'+str(self.m_config['nVar'])+\
                                         '+iVar] '\
                                   '-= kernels::s_m['+str(k)+'];\n')

            # close for loops
            l_sourceFile.write('      }\n'\
                               '    }\n'\
                               '  }\n')

        l_sourceFile.close()






