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

class SpaceTimePredictorGenerator:
    l_config = {}
          
    def __init__(self, i_config):
        l_config = i_config

    
    def generateCode(self):
        self.generatePicardLoop()
        self.generateAverager()
        self.generateExtrapolator()

    
    def generatePicardLoop(self):
        pass

    
    def generateAverager(self):
        pass
    

    def generateExtrapolator(self):
        pass
        # define a sequence of matmul configs