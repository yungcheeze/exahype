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
        # define a sequence of matmul configs
        if(self.l_config['nDim'] == 2):
            #
            # x direction
            #
            l_M = self.l_config['nDof'] * self.l_config['nDof']
            l_N = 2
            l_K = self.l_config['nVar']
            l_lda = l_M
            l_ldb = l_K            
            l_ldc = l_M
            lqhiFCoeff = MatmulConfig(l_M, l_N, l_K, l_lda, l_ldb, l_ldc, 1, 0, 0, 0)
            
            l_commandLineArguments =       lqhiFCoeff.type  + \
                                     ' ' + self.l_pathToLibxsmmGenerator+"/extrapolator.cpp" + \
                                     ' ' + "lqhiFCoeff_" + lqhiFCoeff.baseroutinename + \
                                     ' ' + str(lqhiFCoeff.M) + \
                                     ' ' + str(lqhiFCoeff.N) + \
                                     ' ' + str(lqhiFCoeff.K) + \
                                     ' ' + str(lqhiFCoeff.LDA) + \
                                     ' ' + str(lqhiFCoeff.LDB) + \
                                     ' ' + str(lqhiFCoeff.LDC) + \
                                     ' ' + str(lqhiFCoeff.alpha) + \
                                     ' ' + str(lqhiFCoeff.beta) + \
                                     ' ' + str(lqhiFCoeff.alignment_A) + \
                                     ' ' + str(lqhiFCoeff.alignment_C) + \
                                     ' ' + self.l_config['architecture'] + \
                                     ' ' + "nopf" + \
                                     ' ' + self.l_config['precision'] 
            print(l_commandLineArguments)
            executeLibxsmmGenerator(self.l_pathToLibxsmmGenerator, l_commandLineArguments)                         
            
            #                      
            # y direction
            # 
              
        elif(self.l_config['nDim'] == 3):
            pass
        else:
            print("SpaceTimePredictorGenerator.generateExtrapolator(): nDim not supported")        
            
        
        