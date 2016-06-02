#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Generates the code for the space-time predictor
# for a specific configuration
#--------------------------------------------------------------
#
#
#
import Backend
from MatmulConfig import MatmulConfig
import FunctionSignatures
import Utils


class SpaceTimePredictorGenerator:
    m_config = {}

    # linear/nonlinear
    m_type   = ""

    # for accessing the quadrature weights and so forth
    m_order  = ""

    # total length
    m_luhLength = -1
    m_rhsLength = -1

    # soa size
    m_luhChunkSize = -1

    # padded aos size
    m_structSize = -1


    def __init__(self, i_config, i_numerics):
        self.m_config     = i_config
        self.m_type       = i_numerics
        self.m_order      = str(self.m_config['nDof']-1)
        self.m_structSize = Backend.getSizeWithPadding(self.m_config['nVar'])

        # without padding of lduh, luh
        self.m_luhChunkSize = (self.m_config['nDof']**self.m_config['nDim'])
        self.m_luhLength    = self.m_config['nVar'] * self.m_luhChunkSize
        # with padding of lduh, luh
        #self.m_luhChunkSize    = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        #self.m_luhLength       = self.m_config['nVar'] * self.m_luhChunkSize
        #
        # alternatively
        #self.m_luhChunkSize    = self.m_config['nDof']**self.m_config['nDim']
        #self.m_luhLength       = Backend.getSizeWithPadding(self.m_luhChunkSize)

        self.m_rhsLength = self.m_config['nVar'] * (self.m_config['nDof']**(self.m_config['nDim']+1))


    def __writeHeaderForPicardLoop(self, i_pathToFile):
        l_description = '// Predictor \n'                                                                 \
                        '// K_{1} \\cdot \\hat{q}^{n+1} + K_{\\xi} \\cdot \\hat{F}^{n} = F_0 \\cdot \\hat{u} \n' \
                        '// computed as \n'                                                               \
                        '// \\hat{q}^{n+1} = K_{1}^{-1}[F_0 \\cdot \\hat{u} - K_{\\xi} \\cdot \\hat{F}^{n}] \n'

        l_includeStatement = '#include "string.h"\n'                             \
                             '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/DGMatrices.h"\n'                 \
                             '#include "kernels/GaussLegendreQuadrature.h"\n'    \
                             '#include "kernels/aderdg/optimised/asm_picard.cpp"\n\n'

        l_functionSignature = FunctionSignatures.getPicardLoopSignature(self.m_config['nDim']) + " {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        #l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForPredictor(self, i_pathToFile):
        l_description = '// Compute the time-averaged space-time polynomials (integration in time) \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'              \
                             '#include "kernels/GaussLegendreQuadrature.h"\n'              \
                             '#include "kernels/aderdg/optimised/asm_predictor.cpp"\n\n'

        l_functionSignature = FunctionSignatures.getPredictorSignature()+" {\n"

        l_parameterDocumentation = '// lqh[nVar][nDOFt][nDOFx][nDOFy][nDOFz]      : space-time DOFs\n'                                \
                                   '// lFh[nVar][dim][nDOFx][nDOFy][nDOFz][nDOFt] : space-time DOFs\n'                                \
                                   '// lQhi[nVar][nDOFx][nDOFy][nDOFz]            : time-averaged space-time DOFs\n'                  \
                                   '// lFhi[nVar][nDOF?][nDOF?][nDOF?][dim]       : time-averaged non-linear flux tensor\n'           \
                                   '// where lFhi[nVar][nDOFx][nDOFy][nDOFz][1]\n'                                                    \
                                   '//       lFhi[nVar][nDOFy][nDOFx][nDOFz][2]\n'                                                    \
                                   '//       lFhi[nVar][nDOFz][nDOFy][nDOFx][3]\n\n'

        l_sourceFile = open(i_pathToFile, 'a')
        #l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_parameterDocumentation)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForExtrapolator(self, i_pathToFile):
        l_description = '// Compute the boundary-extrapolated values for Q and F*n \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'      \
                             '#include "kernels/DGMatrices.h"\n'                    \
                             '#include "kernels/aderdg/optimised/asm_extrapolatedPredictor.cpp"\n\n'

        l_functionSignature = FunctionSignatures.getExtrapolatorSignature()+" {\n"

        l_parameterDocumentation = '// lQbnd[nVar][nDOF][nDOF][2*dim]         : boundary-extrapolated data for the state vector\n'    \
                                   '// lFbnd[nVar][nDOFx][nDOFy][nDOFz][2*dim]: the boundary-extrapolated data for the normal flux\n' \
                                   '// lQhi[nVar][nDOFx][nDOFy][nDOFz]        : time-averaged space-time DOFs\n'                      \
                                   '// lFhi[nVar][nDOF?][nDOF?][nDOF?][dim]   : time-averaged non-linear flux tensor\n'               \
                                   '// where lFhi[nVar][nDOFx][nDOFy][nDOFz][1]\n'                                                    \
                                   '//       lFhi[nVar][nDOFy][nDOFx][nDOFz][2]\n'                                                    \
                                   '//       lFhi[nVar][nDOFz][nDOFy][nDOFx][3]\n\n'

        l_sourceFile = open(i_pathToFile, 'a')
        #l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_parameterDocumentation)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForCauchyKovalewski(self, i_pathToFile):
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/DGMatrices.h"\n' \
                             '#include "string.h"\n'                              \
                             '#include "kernels/aderdg/optimised/asm_cauchyKovalewski.cpp"\n\n'
        l_functionSignature = FunctionSignatures.getCauchyKovalewskiSignature()+" {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def generateCode(self):
        if(self.m_type == 'nonlinear'):
            self.__generatePicardLoop()
            self.__generatePredictor()
            self.__generateExtrapolator()
        else:
            self.__generateCauchyKovalewski()
            self.__generateExtrapolator()


    def __generatePicardLoop(self):
        l_filename = "picard.cpp"

        # write #include's and function signature
        self.__writeHeaderForPicardLoop(l_filename)

        l_sourceFile = open(l_filename, 'a')

        # initialisation of lqh and rhs0:
        # (1) lqh(iVar,l,i,j,k) = luh(i,j,k,iVar)
        # (2) rhs0(iVar,i,j,k,l) = weights3(i,j,k) * F0 * luh(iVar,i,j,k)

        nDOF       = str(self.m_config['nDof'])
        blockWidth = str(self.m_config['nDof'] * self.m_structSize)

        # (1) lqh(iVar,l,i,j,k) = luh(i,j,k,iVar);
        if(self.m_config['nDim'] == 2):
            l_sourceFile.write( "  for(int i=0;i<"+nDOF+";i++) {\n"   \
                                "    for(int j=0;j<"+nDOF+";j++) {\n" \
                                "       const int lqh_base_addr = ("+str(self.m_config['nDof']**2)+"*j+"
                                                                    +str(self.m_config['nDof'])   +"*i)*"
                                                                    +str(self.m_structSize)+";\n" \
                                "       const int luh_addr = i + j*"+str(self.m_config['nDof'])+";\n"
                               )
            for iVar in range(0, self.m_config['nVar']):
                l_sourceFile.write("       lqh[lqh_base_addr+"+str(iVar)+"] = luh["+str(iVar*self.m_luhChunkSize)+"+luh_addr];\n")
            l_sourceFile.write( "    }\n"
                                "  }\n"
                              )

        if(self.m_config['nDim'] == 3):
            l_sourceFile.write( "  for(int i=0;i<"+nDOF+";i++) {\n"     \
                                "    for(int j=0;j<"+nDOF+";j++) {\n"   \
                                "      for(int k=0;k<"+nDOF+";k++) {\n" \
                                "         const int lqh_base_addr = ("+str(self.m_config['nDof']**3)+"*k+"
                                                                      +str(self.m_config['nDof']**2)+"*j+"
                                                                      +str(self.m_config['nDof'])   +"*i)*"
                                                                      +str(self.m_structSize)+";\n" \
                                "         const int luh_addr = i + j*"+str(self.m_config['nDof']) +
                                                                "+ k*"+str(self.m_config['nDof']**2)+";\n"
                               )
            for iVar in range(0, self.m_config['nVar']):
                l_sourceFile.write("         lqh[lqh_base_addr+"+str(iVar)+"] = luh["+str(iVar*self.m_luhChunkSize)+"+luh_addr];\n")
            l_sourceFile.write( "      }\n"
                                "    }\n"
                                "  }\n"
                              )



        # 2D/3D
        l_sourceFile.write("  //#pragma omp parallel for\n")
        l_sourceFile.write("  for(int it=0;it<"+str(self.m_config['nDof']**self.m_config['nDim'])+";it++) {\n" \
                           "    const int base_addr = it*"+blockWidth+";\n" \
                           "    for(int l=1;l<"+nDOF+";l++) {\n" \
                           "      memcpy(&lqh[base_addr+l*"+str(self.m_structSize)+"], &lqh[base_addr],"+str(self.m_structSize)+"*sizeof(double));\n"
                           "    }\n" \
                           "  }\n"
                          )


        # (2) rhs0(iVar,i,j,k,l) = weights3(i,j,k) * F0 * luh(iVar,i,j,k)
        # TODO

        # define a sequence of matmul configs
        l_matmulList = []

        # discrete Picard iterations
        # TODO



        # Compute the "derivatives" (contributions of the stiffness matrix)
        # (1) rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) - PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,1,:,j,k,l), Kxi )
        # (2) rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,2,i,:,k,l), Kxi )
        # (3) rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,3,i,j,:,l), Kxi )
        # Incorporate minus sign into dt/dx
        l_sourceFile.write("  // Assume equispaced mesh, dx[0] == dx[1] == dx[2]\n")
        l_sourceFile.write("  const double dtdx = -dt/dx[0];\n\n")

        l_sourceFile.write("  memcpy(rhs,rhs0,"+str(self.m_rhsLength)+"*sizeof(double));\n")

        # (1) rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) + PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,1,:,j,k,l), Kxi )
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    self.m_config['nVar'],                             \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "rhs_x",                                           \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul)

        # write the function call to the driver file
        for i in range(0, self.m_config['nDof']**self.m_config['nDim']):
            l_sourceFile.write(Utils.generateDSCAL("dtdx*kernels::weights3["+str(i)+"]",
                                                   "kernels::Kxi",
                                                   "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
            l_sourceFile.write("  "+l_matmul.baseroutinename
                                   +"(&lFh["+str(i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof'])+"]," \
                                    " &kernels::s_m[0],"  \
                                    " &rhs["+str(i*self.m_config['nVar']*self.m_config['nDof'])+"]);\n\n")



        # (2) rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,2,i,:,k,l), Kxi )


        # (3) rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,3,i,j,:,l), Kxi )


        # lqh(:,:,i,j,k) = MATMUL( rhs(:,i,j,k,:), TRANSPOSE(iK1) )
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']) * self.m_config['nDof']**self.m_config['nDim'], \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    0,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lqh",                                             \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul)

        # write the function call to the driver file
        # note that the DGmatrices.cpp already stores the transpose of iK1
        for i in range(0, self.m_config['nDof']**self.m_config['nDim']):
            l_sourceFile.write(Utils.generateDSCAL("1./kernels::weights3["+str(i)+"]",
                                                   "kernels::iK1",
                                                   "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
            l_sourceFile.write("  "+l_matmul.baseroutinename
                                   +"(&rhs["+str(i*self.m_config['nVar'])+"]," \
                                    " &kernels::s_m[0],"  \
                                    " &lqh["+str(i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof'])+"]);\n\n")

        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)


        # write missing closing bracket
        l_sourceFile.write('}')

        l_sourceFile.close()


    def __generateCauchyKovalewski(self):
        pass


    def __generatePredictor(self):
        l_filename = "predictor.cpp"

        # write #include's and function signature
        self.__writeHeaderForPredictor(l_filename)

        # define a sequence of matmul configs
        l_matmulList = []

        # let's open the file to which we write our function calls to the assembler code
        l_file = open(l_filename, 'a')

        # (1) lqhi(:,i,j,k) = MATMUL(lqh(:,:,i,j,k), wGPN)
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    1,                                                 \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # LDB
                                    self.m_config['nDof'],                             \
                                    # LDC
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    0,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lqhi",                                            \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemv")
        l_matmulList.append(l_matmul)


        # write the function call to the cpp file
        l_iters = self.m_config['nDof'] ** self.m_config['nDim']
        l_baseAddrC = 0
        l_baseAddrA = 0
        for it in range(0, l_iters):
            l_file.write("  "+l_matmul.baseroutinename+'(&lqh['+str(l_baseAddrA)+'],'+  \
                                                        '&kernels::gaussLegendreWeights['+self.m_order+'][0],'+\
                                                        '&lQhi['+str(l_baseAddrC)+'])')
            l_baseAddrC = l_baseAddrC + Backend.getSizeWithPadding(self.m_config['nVar'])
            l_baseAddrA = l_baseAddrA + Backend.getSizeWithPadding(self.m_config['nVar']) * self.m_config['nDof']




        # TODO
        if(self.m_config['nDim'] == 2):
            pass
            # (2) lFhi_x(:,i,j,k) = MATMUL(lFh(:,1,i,j,k,:), wGPN)
            #     lFhi_y(:,j,j,k) = MATMUL(lFh(:,2,i,j,k,:), wGPN)

        elif(self.m_config['nDim'] == 3):
            pass
            # (2) lFhi_x(:,i,j,k) = MATMUL(lFh(:,1,i,j,k,:), wGPN)
            #     lFhi_y(:,j,j,k) = MATMUL(lFh(:,2,i,j,k,:), wGPN)
            #     lFhi_z(:,k,j,i) = MATMUL(lFh(:,3,i,j,k,:), wGPN)

        else:
            print("SpaceTimePredictorGenerator.generatePredictor(): nDim not supported")


        # all matmuls have been collected, now launch code generator backend
        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_file.write('}')
        l_file.close()

    def __generateExtrapolator(self):
        l_filename = "extrapolatedPredictor.cpp"
        
        # write #include's and function signature
        self.__writeHeaderForExtrapolator(l_filename)
                
        # define a sequence of matmul configs        
        l_matmulList = []
               
        # let's open the file to which we write our function calls to the assembler code 
        l_file = open(l_filename, 'a')
       
        if(self.m_config['nDim'] == 2):
            # x-direction
            # (1) lQbnd(:,j,1) = MATMUL(lqhi(:,:,j), FLCoeff)
            # (2) lQbnd(:,j,2) = MATMUL(lqhi(:,:,j), FRCoeff)

            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'],   \
                                     # K                    
                                     self.m_config['nDof'], \
                                     # LDA
                                     1,                     \
                                     # LDB
                                     self.m_config['nDof'], \
                                     # LDC
                                     1,                     \
                                     # alpha 
                                     1,                     \
                                     # beta
                                     0,                     \
                                     # alignment A
                                     0,                     \
                                     # alignment C
                                     0,                     \
                                     # name
                                     "lQbnd",               \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")


            l_matmulList.append(l_matmul)

            # write the function call to the cpp file
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lQhi[0],"                               \
                                                       " &lQbnd[0]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lQhi[0],"                               \
                                                       " &lQbnd["+str(1 * self.m_config['nVar']*self.m_config['nDof'])+"]);\n")

                      
            # y direction
            # (5) lQbnd(:,i,3) = MATMUL(lqhi(:,i,:), FLCoeff)
            # (6) lQbnd(:,i,4) = MATMUL(lqhi(:,i,:), FRCoeff)
            
            # TODO: strided access?
            
            
         
            # x-direction
            # (3) lFbnd(:,j,1) = MATMUL(lFhi_x(:,:,j), FLCoeff)
            # (4) lFbnd(:,j,2) = MATMUL(lFhi_x(:,:,j), FRCoeff)
            #
            # y-direction
            # (7) lFbnd(:,i,3) = MATMUL(lFhi_y(:,:,i), FLCoeff)
            # (8) lFbnd(:,i,4) = MATMUL(lFhi_y(:,:,i), FRCoeff)
            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'],   \
                                     # K                    
                                     self.m_config['nDof'], \
                                     # LDA
                                     1,                     \
                                     # LDB
                                     self.m_config['nDof'], \
                                     # LDC
                                     1,                     \
                                     # alpha 
                                     1,                     \
                                     # beta
                                     0,                     \
                                     # alignment A
                                     0,                     \
                                     # alignment C
                                     0,                     \
                                     # name
                                     "lFbnd",               \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")
            
            l_matmulList.append(l_matmul)
            
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(0 * self.m_config['nVar'] * self.m_config['nDof']**2)+"]," \
                                                       " &lFbnd["+str(0 * self.m_config['nVar']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(0 * self.m_config['nVar'] * self.m_config['nDof']**2)+"]," \
                                                       " &lFbnd["+str(1 * self.m_config['nVar']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(1 * self.m_config['nVar'] * self.m_config['nDof']**2)+"]," \
                                                       " &lFbnd["+str(2 * self.m_config['nVar']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(1 * self.m_config['nVar'] * self.m_config['nDof']**2)+"]," \
                                                       " &lFbnd["+str(3 * self.m_config['nVar']*self.m_config['nDof'])+"]);\n")


        elif(self.m_config['nDim'] == 3):
            # x-direction           
            # (1) lQbnd(:,j,k,1) = MATMUL(lqhi(:,:,j,k), FLCoeff)
            # (2) lQbnd(:,j,k,2) = MATMUL(lqhi(:,:,j,k), FRCoeff)

            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'] * self.m_config['nDof'],   \
                                     # K
                                     self.m_config['nDof'], \
                                     # LDA
                                     1,                     \
                                     # LDB
                                     self.m_config['nDof'], \
                                     # LDC
                                     1,                     \
                                     # alpha
                                     1,                     \
                                     # beta
                                     0,                     \
                                     # alignment A
                                     0,                     \
                                     # alignment C
                                     0,                     \
                                     # name
                                     "lQbnd",               \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")

            l_matmulList.append(l_matmul)
            
            # write the function call to the cpp file
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lQhi[0],"                               \
                                                       " &lQbnd["+str(0 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lQhi[0],"                               \
                                                       " &lQbnd["+str(1 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")



            # y-direction
            # (5) lQbnd(:,i,k,3) = MATMUL(lqhi(:,i,:,k), FLCoeff)
            # (6) lQbnd(:,i,k,4) = MATMUL(lqhi(:,i,:,k), FRCoeff)
            
            # TODO: strided access???
            
            # z-direction
            # (9)  lQbnd(:,i,j,5) = MATMUL(lqhi(:,i,j,:), FLCoeff)
            # (10) lQbnd(:,i,j,6) = MATMUL(lqhi(:,i,j,:), FRCoeff)
            
            # TODO: strided access?
            
            
            

            # x-direction
            # (3) lFbnd(:,j,k,1) = MATMUL(lFhi_x(:,:,j,k), FLCoeff)
            # (4) lFbnd(:,j,k,2) = MATMUL(lFhi_x(:,:,j,k), FRCoeff)
            #
            # y-direction
            # (7) lFbnd(:,i,k,3) = MATMUL(lFhi_y(:,:,i,k), FLCoeff)
            # (8) lFbnd(:,i,k,4) = MATMUL(lFhi_y(:,:,i,k), FRCoeff)
            #
            # z-direction
            # (11) lFbnd(:,i,j,5) = MATMUL(lFhi_z(:,:,i,j), FLCoeff)
            # (12) lFbnd(:,i,j,6) = MATMUL(lFhi_z(:,:,i,j), FRCoeff)
            
            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'] * self.m_config['nDof'],   \
                                     # K
                                     self.m_config['nDof'], \
                                     # LDA
                                     1,                     \
                                     # LDB
                                     self.m_config['nDof'], \
                                     # LDC
                                     1,                     \
                                     # alpha
                                     1,                     \
                                     # beta
                                     0,                     \
                                     # alignment A
                                     0,                     \
                                     # alignment C
                                     0,                     \
                                     # name
                                     "lFbnd",               \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")

            l_matmulList.append(l_matmul)

            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(0 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(0 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(0 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(1 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(1 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(2 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(1 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(3 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FLCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(2 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(4 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::FRCoeff["+self.m_order+"][0]," \
                                                       " &lFhi["+str(2 * self.m_config['nVar'] * self.m_config['nDof']**3)+"],"                               \
                                                       " &lFbnd["+str(5 * self.m_config['nVar']*self.m_config['nDof']*self.m_config['nDof'])+"]);\n")


        else:
            print("SpaceTimePredictorGenerator.generateExtrapolator(): nDim not supported")        


        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_file.write('}')
        l_file.close()