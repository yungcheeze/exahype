#!/bin/env python

from sys import exit

class MatmulConfig:
	# Specification of a dense matrix-matrix multiplication
	#
	#    C       = alpha  *   A   *    B   + beta  *  C
	# (M x N)             (M x K)  (K x N)
	#

	# dense, sparse
	operationType = ''

	baseroutinename = ""
	
	name = ""
	
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
	
	# alignment flags
	alignment_A = 0							# 1 aligned, 0 unaligned  
	alignment_C = 0                         # 1 aligned, 0 unaligned
	
	
	# Constructor
	def __init__(self, M, N, K, LDA, LDB, LDC, alpha, beta, alignment_A, alignment_C, name, operationType='dgemm'):
		if((M > LDC) or (K > LDB) or (M > LDA)):
			print("Incompatible matrix sizes and leading dimensions")
			exit()
		if(alignment_A not in [0,1]):
			print("Something is wrong with the alignment choice of matrix A")
			exit()
		if(alignment_C not in [0,1]):
			print("Something is wrong with the alignment choice of matrix C")
			exit()
			
		self.M = M
		self.N = N
		self.K = K
		self.LDA = LDA
		self.LDB = LDB
		self.LDC = LDC
		self.alpha = alpha
		self.beta = beta
		self.alignment_A = alignment_A
		self.alignment_C = alignment_C
		self.name = name
		self.baseroutinename = operationType+"_"+str(M)+"_"+str(N)+"_"+str(K)
		

	def __repr__(self):
		return "<%s: %s LDA=%s, LDB=%s, LDC=%s, alpha=%s, beta=%s, alignment_A=%s, alignment_C=%s>" \
			 % (self.name, self.baseroutinename, self.LDA, self.LDB, self.LDC, self.alpha, self.beta, self.alignment_A, self.alignment_C)

		