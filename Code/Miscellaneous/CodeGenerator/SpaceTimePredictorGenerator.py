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
from Backend import generateAssemblerCode
from MatmulConfig import MatmulConfig
import FunctionSignatures


class SpaceTimePredictorGenerator:
    m_config = {}

    # for accessing the quadrature weights and so forth
    m_order  = ""
          
    def __init__(self, i_config):
        self.m_config = i_config
        self.m_order  = str(self.m_config['nDof']-1)


    def writeHeaderForPicardLoop(self, i_pathToFile):
        l_description = '// Predictor \n'                                                                 \
                        '// K_{1} \\cdot \\hat{q}^{n+1} + K_{\\xi} \\cdot \\hat{F}^{n} = F_0 \\cdot \\hat{u} \n' \
                        '// computed as \n'                                                               \
                        '// \\hat{q}^{n+1} = K_{1}^{-1}[F_0 \\cdot \\hat{u} - K_{\\xi} \\cdot \\hat{F}^{n}] \n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/DGMatrices.h"\n'                  \
                             '#include "kernels/GaussLegendreQuadrature"\n'       \
                             '#include "kernels/aderdg/optimised/asm_picard.cpp"\n\n'
        
        l_functionSignature = FunctionSignatures.getPicardLoopSignature(self.m_config['nDim']) + " {\n"
        
        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def writeHeaderForPredictor(self, i_pathToFile):
        l_description = '// Compute the time-averaged space-time polynomials (integration in time) \n\n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'           \
                             '#include "kernels/GaussLegendreQuadrature"\n'              \
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
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_parameterDocumentation)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()
                

    def writeHeaderForExtrapolator(self, i_pathToFile):
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
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_parameterDocumentation)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()

    
    def generateCode(self):      
        self.generatePicardLoop()
        self.generatePredictor()
        self.generateExtrapolator()
        
    
    def generatePicardLoop(self):
        pass

    
    def generatePredictor(self):
        l_filename = "predictor.cpp"

        # write #include's and function signature
        self.writeHeaderForPredictor(l_filename)

        # define a sequence of matmul configs
        l_matmulList = []

        # let's open the file to which we write our function calls to the assembler code
        l_file = open(l_filename, 'a')

        if(self.m_config['nDim'] == 2):
            # (1) lqhi(:,i,j,k) = MATMUL(lqh(:,i,j,k,:), wGPN)
            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'] * self.m_config['nDof'], \
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
                                     "lqhi",                \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")

            l_matmulList.append(l_matmul)


            # write the function call to the cpp file
            l_file.write("  //TODO: reorder lqh\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&lQhi[0], "+\
                                                       "&kernels::gaussLegendreWeights["+self.m_order+"][0],"\
                                                       "&lqh[0]);\n")
   


            # (2) lFhi_x(:,i,j,k) = MATMUL(lFh(:,1,i,j,k,:), wGPN)
            #     lFhi_y(:,j,j,k) = MATMUL(lFh(:,2,i,j,k,:), wGPN)
            #     lFhi_z(:,k,j,i) = MATMUL(lFh(:,3,i,j,k,:), wGPN)
            # TODO


        elif(self.m_config['nDim'] == 3):
            # (1) lqhi(:,i,j,k) = MATMUL(lqh(:,i,j,k,:), wGPN)
            l_matmul = MatmulConfig( # M
                                     1,                     \
                                     # N
                                     self.m_config['nVar'] * self.m_config['nDof'] * self.m_config['nDof'] * self.m_config['nDof'], \
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
                                     "lqhi",                \
                                     # prefetching
                                     "nopf",                \
                                     # type
                                     "gemv")

            l_matmulList.append(l_matmul)

            # write the function call to the cpp file
            l_file.write("  //TODO: reorder lqh\n")
            l_file.write("  "+l_matmul.baseroutinename+"(&kernels::gaussLegendreWeights["+self.m_order+"][0]," \
                                                       " &lqh[0],"                                             \
                                                       " &lQhi[0]\n")


            # (2) lFhi(:,iDim,i,j,k) = MATMUL(lFh(:,iDim,i,j,k,:), wGPN)
            # TODO

        else:
            print("SpaceTimePredictorGenerator.generatePredictor(): nDim not supported")


        # all matmuls have been collected, now launch code generator backend
        generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_file.write('}')
        l_file.close()

    def generateExtrapolator(self):
        l_filename = "extrapolatedPredictor.cpp"
        
        # write #include's and function signature
        self.writeHeaderForExtrapolator(l_filename)
                
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


        generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_file.write('}')
        l_file.close()        
        