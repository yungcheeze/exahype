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
from Backend import executeLibxsmmGenerator
from MatmulConfig import MatmulConfig


class SpaceTimePredictorGenerator:
    l_config = {}
    l_pathToLibxsmmGenerator = ""
          
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
        
        if(self.l_config['nDim']==2):                        
            l_functionSignature = "template<void PDEFlux2d(const double * const Q, double * f, double * g)>\n" \
                                  "void kernels::aderdg::optimised::picardLoop( \n"                            \
                                  "  double * lQi, \n"                                                         \
                                  "  double * lFi, \n"                                                         \
                                  "  const double * const luh, \n"                                             \
                                  "  const tarch::la::Vector<DIMENSIONS,double> &dx\n"                         \
                                  "  const double dt \n"                                                       \
                                  ") {\n"
        elif(self.l_config['nDim']==3):
            l_functionSignature = "template<void PDEFlux3d(const double * const Q, double * f, double * g, double * h)>\n" \
                                  "void kernels::aderdg::optimised::picardLoop( \n"                                        \
                                  "  double * lQi, \n"                                                                     \
                                  "  double * lFi, \n"                                                                     \
                                  "  const double * const luh, \n"                                                         \
                                  "  const tarch::la::Vector<DIMENSIONS,double> &dx\n"                                     \
                                  "  const double dt \n"                                                                   \
                                  ") {\n"
        else:
            print("SpaceTimePredictorGenerator.writeHeaderForPicardLoop(): nDim not supported") 
        
        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def writeHeaderForPredictor(self, i_pathToFile):
        l_description = '// Compute the time-averaged space-time polynomials (integration in time) \n\n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/GaussLegendreQuadrature"\n\n'     
                                
        l_functionSignature = "void kernels::aderdg::optimised::predictor( \n"    \
                              "  double * lQhi, \n"                               \
                              "  double * lFhi, \n"                               \
                              "  double * lQhi, \n"                               \
                              "  double * lFi \n"                                 \
                              ") {\n"
        
        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()
                

    def writeHeaderForExtrapolator(self, i_pathToFile):
        l_description = '// Compute the boundary-extrapolated values for Q and F*n \n\n'
        
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/DGMatrices.h"\n\n'
                                
        l_functionSignature = "void kernels::aderdg::optimised::extrapolator( \n" \
                              "  double * lQbnd, \n"                              \
                              "  double * lFbnd, \n"                              \
                              "  double * lqhi, \n"                               \
                              "  double * lFhi \n"                                \
                              ") {\n"
                              
        l_parameterDocumentation = '// lQbnd[nVar][nDOF][nDOF][nFace]       : boundary-extrapolated data for the state vector\n' \
                                   '// lFbnd[nVar][nDOFx][nDOFy][nDOFz][dim]: the boundary-extrapolated data for the normal flux\n' \
                                   '// lqhi[nVar][nDOFx][nDOFy][nDOFz]      : time-averaged space-time DOFs\n'\
                                   '// lFhi[nVar][nDOFx][nDOFy][nDOFz][dim] : time-averaged non-linear flux tensor\n\n'      
        
        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_parameterDocumentation)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()

    
    def generateCode(self, i_pathToLibxsmmGenerator):
        self.l_pathToLibxsmmGenerator = i_pathToLibxsmmGenerator
        
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
        
        # we do classical DGEMM: C = A * B
        l_alpha = 1
        l_beta  = 0
        
        # let's open the file to which we write our function calls to the assembler code 
        l_file = open(l_filename, 'a')
       
        if(self.l_config['nDim'] == 2):       
            #
            # x-direction
            #
            
            # lqhi * FLCoeff, unpadded
            # LDA, LDB, LDC, alpha, beta, alignment_A, alignment_C, baseroutinename
            l_matmulList.append(MatmulConfig(self.l_config['nVar'],1, self.l_config['nDof'], \
                                             self.l_config['nVar'], self.l_config['nDof'], self.l_config['nVar'], \
                                             l_alpha, l_beta, 0, 0, "lQbnd"))    
            
            # number of entries between two flux matrices, or, equivalently, the number of face DOFs
            l_offset = self.l_config['nVar']*self.l_config['nDof']
            
            # write the function calls to the cpp file
            for l_matmul in l_matmulList:
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
                
                                  
            #                      
            # y direction
            # 
            
            
            """
            l_commandLineArguments =       lqhiFLCoeff.type  + \
                                     ' ' + self.l_pathToLibxsmmGenerator+"/extrapolatedPredictor.cpp" + \
                                     ' ' + lqhiFLCoeff.baseroutinename + \
                                     ' ' + str(lqhiFLCoeff.M) + \
                                     ' ' + str(lqhiFLCoeff.N) + \
                                     ' ' + str(lqhiFLCoeff.K) + \
                                     ' ' + str(lqhiFLCoeff.LDA) + \
                                     ' ' + str(lqhiFLCoeff.LDB) + \
                                     ' ' + str(lqhiFLCoeff.LDC) + \
                                     ' ' + str(lqhiFLCoeff.alpha) + \
                                     ' ' + str(lqhiFLCoeff.beta) + \
                                     ' ' + str(lqhiFLCoeff.alignment_A) + \
                                     ' ' + str(lqhiFLCoeff.alignment_C) + \
                                     ' ' + self.l_config['architecture'] + \
                                     ' ' + "nopf" + \
                                     ' ' + self.l_config['precision'] 
            print(l_commandLineArguments)
            executeLibxsmmGenerator(self.l_pathToLibxsmmGenerator, l_commandLineArguments)
            """   
            
              
        elif(self.l_config['nDim'] == 3):
            pass
        else:
            print("SpaceTimePredictorGenerator.generateExtrapolator(): nDim not supported")        
            
        # write missing closing bracket    
        l_file.write('}')
        l_file.close()        
        