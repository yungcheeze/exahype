SUBROUTINE BoundaryConditions 
    USE typesDef
    ! Local variables 
    INTEGER :: iFace 
    REAL    :: j,k
    REAL    :: Qbc(nVar),Fbc(nVar,d),Vbc(nVar)  
    INTEGER :: ibc(1)  
    !
#ifdef ELASTICITY
    Qbc = 0. 
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use either simple extrapolation (copy from inside the domain) 
        ! or impose a constant value 
        !qbc = (/ 1., 0., 0., 0., 2.5 /) 
        IF(Face(iFace)%Left.EQ.0) THEN
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qL(:,j,k) = Face(iFace)%qR(:,j,k)
                    IF( ibc(1)==2 ) THEN
                        Face(iFace)%qL( (/2,4,5/),j,k) = -Face(iFace)%qL((/2,4,5/),j,k)
                    !    Face(iFace)%qL(:,j,k) = qbc
                    ELSE 
                        Face(iFace)%qL(:,j,k) = qbc 
                    ENDIF 
                    !Face(iFace)%FL(:,j,k) = Face(iFace)%FR(:,j,k)
                    !Face(iFace)%qL(:,j,k) = qbc  
                    !CALL PDEFlux(Fbc,Face(iFace)%qL(:,j,k),Face(iFace)%paramL(:,j,k)) 
                    !Face(iFace)%FL(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    Face(iFace)%qR(:,j,k) = Face(iFace)%qL(:,j,k) 
                    IF( ibc(1)==2 ) THEN
                        Face(iFace)%qR( (/2,4,5/),j,k) = -Face(iFace)%qR((/2,4,5/),j,k)
                    !    Face(iFace)%qR(:,j,k) = qbc
                    ELSE 
                        Face(iFace)%qR(:,j,k) = qbc 
                    ENDIF 
                    !Face(iFace)%FR(:,j,k) = Face(iFace)%FL(:,j,k)
                    !Face(iFace)%qR(:,j,k) = qbc  
                    !CALL PDEFlux(Fbc,Face(iFace)%qR(:,j,k),Face(iFace)%paramR(:,j,k)) 
                    !Face(iFace)%FR(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO    
#else
    ! Fix boundary data  
    Vbc = (/ 1., 1., 1., 0., 1. /)    ! primitive variables     
    CALL PDEPrim2Cons(qBC,Vbc)        ! convert into conservative variables    
    DO iFace = 1, nFace
        ! Here, we need to take care of the boundary conditions 
        ! For the moment, we use either simple extrapolation (copy from inside the domain) 
        ! or impose a constant value 
        IF(Face(iFace)%Left.EQ.0) THEN
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    !Face(iFace)%qL(:,j,k) = Face(iFace)%qR(:,j,k)
                    !Face(iFace)%FL(:,j,k) = Face(iFace)%FR(:,j,k)
                    Face(iFace)%qL(:,j,k) = qbc 
                    CALL PDEFlux(Fbc,Face(iFace)%qL(:,j,k),Face(iFace)%paramL(:,j,k)) 
                    Face(iFace)%FL(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF
        IF(Face(iFace)%Right.EQ.0) THEN             
            ibc = MAXLOC( ABS(Face(iFace)%nv(:)) )
            ibc = SIGN(1.0, Face(iFace)%nv(ibc(1)) ) * ibc  
            DO k = 1, nDOF(3)
                DO j = 1, nDOF(2) 
                    !Face(iFace)%qR(:,j,k) = Face(iFace)%qL(:,j,k) 
                    !Face(iFace)%FR(:,j,k) = Face(iFace)%FL(:,j,k)
                    Face(iFace)%qR(:,j,k) = qbc 
                    CALL PDEFlux(Fbc,Face(iFace)%qR(:,j,k),Face(iFace)%paramR(:,j,k)) 
                    Face(iFace)%FR(:,j,k) = MATMUL(Fbc(:,:), Face(iFace)%nv) 
                ENDDO
            ENDDO 
        ENDIF            
    ENDDO    

#endif 
    ! 
END SUBROUTINE BoundaryConditions 
    
    