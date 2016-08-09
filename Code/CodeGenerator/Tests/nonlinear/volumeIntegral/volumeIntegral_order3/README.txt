Test parameters
***************

Kernel tested: Volume Integral (2d)
Order: 3
Non-Linear
nVar = 5
Default padding = 3


Test input
**********

* Fhi(nVar,d,nDOF(1),nDOF(2),nDOF(3))
	with
		nVar = 5
		d = 2
		nDOF(1) = nDOF(2) = 4
		nDOF(3) = 1
	Input format:
		for(d) 
			for(nVar)
                for(nD1) 
                  Fhi(nVar,d,nD1,:,1)
				
