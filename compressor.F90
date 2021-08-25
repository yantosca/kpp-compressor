program main

  USE GCKPP_MODEL
  USE GCKPP_INITIALIZE
  
  IMPLICIT NONE

  INTEGER                :: ICNTRL(20), IERR, I, N
  INTEGER                :: ISTATUS(20)
  REAL(dp)               :: RCNTRL(20)
  REAL(dp)               :: RSTATE(20)
  REAL(dp)               :: T, TIN, TOUT, start, end
  REAL :: sumtime

  ! COMPACTION
  INTEGER, ALLOCATABLE :: REMOVE(:)    ! Vector of species indexes to be removed from full mechanism
!  INTEGER, ALLOCATABLE :: cLU_IROW(:)  ! Compacted ROW indexes
!  INTEGER, ALLOCATABLE :: cLU_ICOL(:)  ! Compacted COL indexes
!  INTEGER, ALLOCATABLE :: cLU_CROW(:)  ! Compacted compressed row vector
!  INTEGER, ALLOCATABLE :: cLU_DIAG(:)  ! Compacted DIAG indexes
  INTEGER :: NRMV                      ! Number of species removed (size(REMOVE))
!  INTEGER :: cNONZERO                  ! Number of non-zero entries in compacted Jacobian
!  INTEGER :: rNVAR                     ! Reduced number of variable species (the name "cNVAR" is already taken by KPP)
  
  INTEGER :: idx, S

  ALLOCATE(DO_SLV(NVAR+1))     ! "VAR"iable"L"ogical, used to control KppSolve()
  ALLOCATE(DO_JVS(NONZERO)) ! Yes/No, compute term, used to control Jac_SP()

  ! Temporoary settings.
  DO_SLV     = 1 ! Compute all in KppSolve()
  DO_JVS = 1 ! Compute all in JacSP()
  cNONZERO = 0 ! Initialize

  ! -------------------------------------------------------------------------- !
  ! 1. Run the full mechanism

  call fullmech()

  ! -------------------------------------------------------------------------- !
  ! 2. Reconstruct the sparse data for a reduced mechanism
  ! e.g. compact the Jacobian

  NRMV     = 1 ! Remove 1 species for testing. This would be determined online.
  rNVAR    = NVAR-NRMV ! Number of active species in the reduced mechanism

  ! ALLOCATE
  ALLOCATE(REMOVE(NRMV))

  ! -- remove row & column
  ! -- -- Which species are zeroed?
  REMOVE(1) = ind_B ! Species B

  ! -- DO_SLV and DO_JVS will not change size (remain NVAR & NONZERO)
  !    But the appropriate elements are set to zero, so the appropriate terms
  !    in KppSolve() and Jac_SP() are not computed

  DO_SLV(REMOVE(:)) = 0  

  ! -- -- Loop through vectors
  ! -- -- -- Count the number of nonzero elements in the reduced Jacobian
  DO S = 1,NRMV
  DO i = 1,NONZERO !NONZERO is the size of LU_ICOL & LU_IROW
     IF (LU_IROW(i).ne.REMOVE(S).and.LU_ICOL(i).ne.REMOVE(S)) cNONZERO = cNONZERO+1
  ENDDO
  ENDDO

  ! -- -- -- Allocate the new Jacobian elements based on new non-zero elements
  ALLOCATE(cLU_IROW(cNONZERO))
  ALLOCATE(cLU_ICOL(cNONZERO))
  ALLOCATE(JVS_MAP(cNONZERO))
  ALLOCATE(cLU_CROW(rNVAR+1))
  ALLOCATE(cLU_DIAG(rNVAR+1))

  ! -- -- -- recompute LU_IROW
  ! -- -- -- recompute LU_ICOL
  DO S = 1,NRMV
  idx = 0
  DO i = 1,NONZERO !NONZERO is the size of LU_ICOL & LU_IROW
     ! Remove row & column elements of removed species
     ! Populate cLU_IROW & cLU_ICOL
     IF (LU_IROW(i).ne.REMOVE(S).and.LU_ICOL(i).ne.REMOVE(S)) THEN
        idx=idx+1
        IF (LU_IROW(i).gt.S) THEN
           cLU_IROW(idx) = LU_IROW(i)-S
           cLU_ICOL(idx) = LU_ICOL(i)-S
           JVS_MAP(idx)  = i
        ELSE
           cLU_IROW(idx) = LU_IROW(i)
           cLU_ICOL(idx) = LU_ICOL(i)
           JVS_MAP(idx)  = i
        ENDIF
     ENDIF
     ! Toggle DO_JVS()
     IF (LU_IROW(i).eq.REMOVE(S).or.LU_ICOL(i).eq.REMOVE(S)) THEN
        DO_JVS(i) = 0 ! Turn off this term in Jac_SP()
     ENDIF
  ENDDO
  ENDDO

  cLU_CROW(1) = 1 ! 1st index = 1
  cLU_DIAG(1) = 1 ! 1st index = 1

  ! -- -- -- recompute LU_CROW
  S = 2
  DO i = 2,cNONZERO
     IF (cLU_IROW(i).ne.cLU_IROW(i-1).and.S.le.rNVAR) THEN
        cLU_CROW(S) = i
        S=S+1
     ENDIF
  ENDDO

  ! -- -- -- recompute LU_DIAG
  S = 1
  DO i = 2,cNONZERO
     IF (cLU_IROW(i).eq.cLU_ICOL(i)) THEN
        S=S+1
        cLU_DIAG(S) = i 
     ENDIF
  ENDDO

  cLU_CROW(rNVAR+1) = cNONZERO+1
  cLU_DIAG(rNVAR+1) = cLU_DIAG(rNVAR)+1

!  write(*,*) 'cNONZERO: ', cNONZERO
  write(*,*) ' '
  write(*,*) 'Species Info:'
  write(*,*) '---------------'
  do i=1,nvar
     write(*,'(a,i3)') " Species "//trim(spc_names(i))//" has index: ", i
  enddo
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'Full-mech sparse data: '
  write(*,*) 'LU_IROW: ', LU_IROW
  write(*,*) 'LU_ICOL: ', LU_ICOL
  write(*,*) 'LU_CROW: ', LU_CROW
  write(*,*) 'LU_DIAG: ', LU_DIAG
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'Compacted sparse data: '
  write(*,*) '-- removes species ' // SPC_NAMES(REMOVE(:))! // ' with index ', REMOVE(:)
  write(*,*) 'cLU_IROW: ', cLU_IROW
  write(*,*) 'cLU_ICOL: ', cLU_ICOL
  write(*,*) 'cLU_CROW: ', cLU_CROW
  write(*,*) 'cLU_DIAG: ', cLU_DIAG
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'JVS_MAP ensures that the right JVS values are indexed in KppDecomp()'
  write(*,*) 'JVS_MAP: ', JVS_MAP
  write(*,*) ' '
  write(*,*) 'DO_SLV controls the terms that will be computed in KppSolve(): 1=compute, 0=skip'
  write(*,*) 'DO_SLV:  ', DO_SLV
  write(*,*) ' '
  write(*,*) 'DO_JVS controls the terms that will be computed in Jac_SP(): 1=compute, 0=skip'
  write(*,*) 'DO_JVS:  ', DO_JVS

  ! -------------------------------------------------------------------------- !
  ! 3. Run the compacted mechanism
  ! - Need to somehow pass the compacted vectors to KPP
  ! - Should be pretty straight-forward
  ! - Will need to respond to new 'parameters' cNVAR, cNONZERO
  ! - the compacted mechanism will still process the full species vector, 
  !   C(NVAR), not c(rNVAR)

  ! with DO_JVS() and DO_SLV set above, this will 'try' to run a reduced mech, but will
  ! dump errors since it still runs the full mech. Goal is to make it properly run the 
  ! reduced mech.
  ! call compactedmech()

  ! -------------------------------------------------------------------------- !
  ! 4. (optional) Calculate the error norm per Santillana et al. (2010) and
  !    Shen et al. (2020)
  !    --- this is valuable if the species removed are selected per
  !        a reduction criterion. Requires a properly posed mechanism.

  ! -------------------------------------------------------------------------- !

CONTAINS

  subroutine fullmech()
    IMPLICIT NONE
    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0.0_dp            ! Rosenbrock output 
    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RSTATE    = 0.0_dp            ! Rosenbrock output
    ICNTRL    = 0
    ICNTRL(1) = 1
    ICNTRL(2) = 0	
    ICNTRL(3) = 4
    ICNTRL(7) = 1
    
    ! Tolerances
    ATOL      = 1e-2_dp    
    RTOL      = 1e-2_dp
    
    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 3600._dp
    TEMP = 298.
    
    write(*,*) ' '
    write(*,*) 'Running the full mechanism'
    
    call cpu_time(start)
    sumtime   = 0.
    ! Initialize
    call Initialize()
    
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,20
       DO N=1,1000
          ! Set C
          C(1:NVAR)   = 1e8
          VAR(1:NVAR) = C(1:NVAR)
          
          ! Set RCONST
          R = 1.
          !        R(2) = 1.
          call Update_RCONST()
          
          ! Integrate
          CALL Integrate( TIN,    TOUT,    ICNTRL,      &
               RCNTRL, ISTATUS, RSTATE, IERR )
          C(1:NVAR)       = VAR(:)
          
       ENDDO
       call cpu_time(end)
       sumtime = sumtime+end-start
    ENDDO
    write(*,*) "Average integration time: ", sumtime/20.
    !  write(*,*) ' '
    write(*,*) '---------------'
  end subroutine fullmech

  subroutine compactedmech()
    IMPLICIT NONE
    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0.0_dp            ! Rosenbrock output 
    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RSTATE    = 0.0_dp            ! Rosenbrock output
    ICNTRL    = 0
    ICNTRL(1) = 1
    ICNTRL(2) = 0	
    ICNTRL(3) = 4
    ICNTRL(7) = 1
    
    ! Tolerances
    ATOL      = 1e-2_dp    
    RTOL      = 1e-2_dp
    
    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 3600._dp
    TEMP = 298.
    
    write(*,*) ' '
    write(*,*) 'Running the reduced mechanism'
    
    call cpu_time(start)
    sumtime   = 0.
    ! Initialize
    call Initialize()
    
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,20 ! Iterate to generate average comp time
       DO N=1,1000
          ! Set C
          C(1:NVAR)   = 1e8
          VAR(1:NVAR) = C(1:NVAR)
          
          ! Set RCONST
          R = 1.
          !        R(2) = 1.
          call Update_RCONST()
          
          ! Integrate
          CALL Integrate( TIN,    TOUT,    ICNTRL,      &
               RCNTRL, ISTATUS, RSTATE, IERR )
          C(1:NVAR)       = VAR(:)
          
       ENDDO
       call cpu_time(end)
       sumtime = sumtime+end-start
    ENDDO
    write(*,*) "Average integration time: ", sumtime/20.
    write(*,*) '---------------'
  end subroutine compactedmech

end program main
