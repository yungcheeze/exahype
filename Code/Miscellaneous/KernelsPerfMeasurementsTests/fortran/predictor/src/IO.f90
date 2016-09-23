SUBROUTINE WriteData 
  USE typesDef   
  USE ISO_C_BINDING
  IMPLICIT NONE 
  !include 'tecio.f90' 
  CHARACTER(LEN=200) :: Filename,Title,ScratchDir, VarString   
  CHARACTER(LEN=10)  :: Name(nVar) 
  CHARACTER(LEN=10)  :: cmyrank 
  INTEGER            :: i,j,ii,jj,c,nc,iRet,iDim,iErr 
  REAL               :: QN(nVar),VN(nVar),Vav(nVar),Vmin(nVar),ldx(d),lx0(d),lxb(d),ux(nVar),uy(nVar)
  REAL               :: LocNode(nVar,(N+1)**nDim), GradNode(nVar,(N+1)**nDim,d),  xvec(d) 
  INTEGER            :: nSubNodes, nPlotElem, nSubPlotElem, nRealNodes, ZoneType, nVertex, nDOFs   
  REAL(8)            :: loctime 
  INTEGER*4, POINTER :: NData(:,:)  
  REAL*4, POINTER    :: DataArray(:,:),TempArray(:) 
  INTEGER*4          :: visdouble
  REAL*4             :: Test 
  POINTER   (NullPtr,Null)
  Integer*4 Null(*)
  !
  visdouble = 0
  nVertex = 2**nDim 
  WRITE(FileName,'(a,a1,i8.8,a)') TRIM(BaseFile),'-', timestep, '.plt'
  !PRINT *, ' Writing data to file ', TRIM(FileName) 
  !
  NullPtr = 0 
  nPlotElem =  0
  nSubPlotElem = 0
  nRealNodes = 0  
  DO i = 1, nElem
     nPlotElem = nPlotElem + 1
     nSubPlotElem = nSubPlotElem + N**nDim  
     nSubNodes = (N+1)**nDim  
     nRealNodes = nRealNodes + nSubNodes 
  ENDDO
  ALLOCATE(NData(nVertex,nSubPlotElem))  
  ALLOCATE(DataArray(nRealNodes,nDim+nVar+1))  
  c  = 0 
  nc = 0 
  DO i = 1, nElem
    DO j = 1, N**nDim 
        c = c + 1  
        NData(:,c) = nc + subtri(1:nVertex,j)
    ENDDO 
    nc = nc + (N+1)**nDim   
  ENDDO
  nDOFs = PRODUCT(nDOF(1:nDim)) 
  c = 0 
  DO i = 1, nElem      
    LocNode = MATMUL( RESHAPE( uh(:,:,:,:,i), (/ nVar, nDOFs /) ), SubOutputMatrix(1:nDOFs,1:(N+1)**nDim) ) 
    lx0 = x(:,tri(1,i)) 
    DO j = 1, (N+1)**nDim  
        QN(:) = LocNode(:,j) 
        xvec = lx0 + allsubxi(:,j)*dx 
        CALL PDECons2Prim(VN,QN,iErr)
        c = c + 1 
        DataArray(c,:) = (/ xvec(1:nDim), VN, REAL(i) /)   
    ENDDO
  ENDDO

  WRITE(Title,'(a,f9.4,a)') 'Time t = ', time, ''//C_NULL_CHAR  
  WRITE(ScratchDir,'(a)') '.'//C_NULL_CHAR 
  SELECT CASE(nDim)
  CASE(1)
   WRITE(VarString,*) 'x ' 
   ZoneType = 1 ! FEM Line seg   
  CASE(2)
   WRITE(VarString,*) 'x y ' 
   ZoneType = 3 ! FEM Quad  
  CASE(3)
   WRITE(VarString,*) 'x y z ' 
   ZoneType = 5 ! FEM Brick 
  END SELECT 
  CALL PDEVarName(Name)  
  DO i = 1, nVar        
     WRITE(VarString,'(a,a,a,a)') TRIM(VarString), ' ', TRIM(Name(i)) , ' '   
  ENDDO
  WRITE(VarString,'(a,a)') TRIM(VarString), ' iE ' 

  loctime = time 


  DEALLOCATE(NData,DataArray)  

END SUBROUTINE WriteData 
    

    