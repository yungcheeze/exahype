Test parameters
***************

Kernel tested: Volume Integral (2d)
Order: 3
Non-Linear
nVar = 5
Default padding = 3


Test input
**********

* lduh(nVar,nDOF(1),nDOF(2),nDOF(3))
	with
		nVar = 5
		d = 2
		nDOF(1) = nDOF(2) = 4
		nDOF(3) = 1
	Input format:
		DO k=1,nDOF(3)
          DO j=1,nDOF(2)
            DO i=1,nDOF(1)
              READ(99,*) lduh(:,i,j,k)
            ENDDO
          ENDDO
        ENDDO
        
* lFbnd(nVar,6, nDOF(2),nDOF(3))
	with
		nVar = 5
		d = 2
		nDOF(1) = nDOF(2) = 4
		nDOF(3) = 1
		6 = nFace (last 2 unused in 2d)
	Input format:
		DO k=1,nDOF(3)
		  DO j=1,nDOF(2)
			DO iVar=1,nVar
			  READ(99,*) lFbnd(iVar,:,j,k)
			ENDDO
		  ENDDO
		ENDDO