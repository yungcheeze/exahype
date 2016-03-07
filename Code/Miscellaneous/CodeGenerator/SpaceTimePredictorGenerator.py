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
    l_config = {}
          
    def __init__(self, i_config):
        self.l_config = i_config


    def writeHeaderForPicardLoop(self, i_pathToFile):
        l_description = '// Predictor \n'                                                                 \
                        '// K_{1} \\cdot \\hat{q}^{n+1} + K_{\\xi} \\cdot \\hat{F}^{n} = F_0 \\cdot \\hat{u} \n' \
                        '// computed as \n'                                                               \
                        '// \\hat{q}^{n+1} = K_{1}^{-1}[F_0 \\cdot \\hat{u} - K_{\\xi} \\cdot \\hat{F}^{n}] \n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/DGMatrices.h"\n'                  \
                             '#include "kernels/GaussLegendreQuadrature"\n\n'     
        
        l_functionSignature = FunctionSignatures.getPicardLoopSignature(self.l_config['nDim']) + " {\n"
        
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
        
        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()
                

    def writeHeaderForExtrapolator(self, i_pathToFile):
        l_description = '// Compute the boundary-extrapolated values for Q and F*n \n\n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'      \
                             '#include "kernels/DGMatrices.h"\n'                    \
                             '#include "kernels/aderdg/optimised/asm_extrapolatedPredictor.cpp"\n\n'
                                
        l_functionSignature = FunctionSignatures.getExtrapolatorSignature()+" {\n"
                              
        l_parameterDocumentation = '// lQbnd[nVar][nDOF][nDOF][nFace]       : boundary-extrapolated data for the state vector\n'    \
                                   '// lFbnd[nVar][nDOFx][nDOFy][nDOFz][dim]: the boundary-extrapolated data for the normal flux\n' \
                                   '// lqhi[nVar][nDOFx][nDOFy][nDOFz]      : time-averaged space-time DOFs\n'                      \
                                   '// lFhi[nVar][nDOF?][nDOF?][nDOF?][dim] : time-averaged non-linear flux tensor\n'               \
                                   '// where lFhi[nVar][nDOFx][nDOFy][nDOFz][1]\n'                                                  \
                                   '//       lFhi[nVar][nDOFy][nDOFx][nDOFz][2]\n'                                                  \
                                   '//       lFhi[nVar][nDOFy][nDOFx][nDOFz][3]\n\n'                                                  
        
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
        pass
    

    def generateExtrapolator(self):
        l_filename = "extrapolatedPredictor.cpp"
        
        # write #include's and function signature
        self.writeHeaderForExtrapolator(l_filename)
                
        # define a sequence of matmul configs        
        l_matmulList = []        
               
        # let's open the file to which we write our function calls to the assembler code 
        l_file = open(l_filename, 'a')
       
        if(self.l_config['nDim'] == 2):
            # ---------------------------------
            # unpadded, dgemv variant
            # ---------------------------------
            
                   
            # 
            # x-direction
            # 
            
            # we do classical DGEMM: C = A * B
            l_matmul = MatmulConfig( # M
                                     self.l_config['nVar'], \
                                     # N
                                     1,                     \
                                     # K                    
                                     self.l_config['nDof'], \
                                     # LDA
                                     self.l_config['nVar'], \
                                     # LDB
                                     self.l_config['nDof'], \
                                     # LDC
                                     self.l_config['nVar'], \
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
                                     "nopf")
            
                     
            # number of entries between two flux matrices, or, equivalently, the number of face DOFs
            l_offset = self.l_config['nVar']*self.l_config['nDof']
            
            # write the function calls to the cpp file
            for j in range(self.l_config['nDof']):
                # lQbnd(:,1,j,k) = MATMUL( lqhi(:,:,j,k),   FLCoeff )   ! left
                l_file.write("  "+l_matmul.baseroutinename+"(&lqhi["+str(j*self.l_config['nVar']*self.l_config['nDof'])+"], "+\
                                                           "FLCoeff,"\
                                                           " &lQbnd["+str(j*self.l_config['nVar'])+"]);\n")
                # lQbnd(:,2,j,k) = MATMUL( lqhi(:,:,j,k),   FRCoeff )   ! right
                l_file.write("  "+l_matmul.baseroutinename+"(&lqhi["+str(j*self.l_config['nVar']*self.l_config['nDof'])+"], "+\
                                                           "FRCoeff,"\
                                                           " &lQbnd["+str(j*self.l_config['nVar']+l_offset)+"]);\n")
                # lFbnd(:,1,j,k) = MATMUL( lFhi(:,1,:,j,k), FLCoeff )   ! left
                l_file.write("  "+l_matmul.baseroutinename+"(&lFhi["+str(j*self.l_config['nVar']*self.l_config['nDof'])+"], "+\
                                                           "FLCoeff,"\
                                                           " &lFbnd["+str(j*self.l_config['nVar'])+"]);\n")
                # lFbnd(:,2,j,k) = MATMUL( lFhi(:,1,:,j,k), FRCoeff )   ! right
                l_file.write("  "+l_matmul.baseroutinename+"(&lFhi["+str(j*self.l_config['nVar']*self.l_config['nDof'])+"], "+\
                                                           "FRCoeff,"\
                                                           " &lFbnd["+str(j*self.l_config['nVar']+l_offset)+"]);\n")
                
            l_matmulList.append(l_matmul)
                                                      
            #                      
            # y direction
            # 
            
            generateAssemblerCode("asm_"+l_filename, l_matmulList)
            
              
        elif(self.l_config['nDim'] == 3):
            pass
        else:
            print("SpaceTimePredictorGenerator.generateExtrapolator(): nDim not supported")        
            
        # write missing closing bracket    
        l_file.write('}')
        l_file.close()        
        