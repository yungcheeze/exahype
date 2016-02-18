#!/bin/env python

from sys import exit

class MatmulConfig:
	# Specification of a dense matrix-matrix multiplication
	#
	#    C       = alpha  *   A   *    B   + beta  *  C
	# (M x N)             (M x K)  (K x N)
	#

	baseroutinename = ""
	
	# dimension of the matrices
	M = -1
	N = -1
	K = -1
	
	# leading dimension of A, B, and C 
	LDA = -1
	LDB = -1
	LDC = -1
	
	# scalars
	alpha = 1;								# -1, 1
	beta  = 1;								#  0, 1
	
	
	
	# Constructor
	def __init__(self, M, N, K, LDA, LDB, LDC, alpha, beta):
		if((M > LDC) or (K > LDB) or (M > LDA)):
			print("Incompatible matrix sizes and leading dimensions")
			exit()
			
		self.M   	= M
		self.N   	= N
		self.K   	= K
		self.LDA 	= LDA
		self.LDB 	= LDB
		self.LDC 	= LDC
		self.alpha 	= alpha
		self.beta 	= beta
		self.baseroutinename = "dense_"+str(M)+"_"+str(N)+"_"+str(K)
		

		