program main

  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE SETQUANTS

  IMPLICIT NONE

  INTEGER                :: ICNTRL(20), IERR, I, N
  INTEGER                :: ISTATUS(20)
  REAL(dp)               :: RCNTRL(20)
  REAL(dp)               :: RSTATE(20)
  REAL(dp)               :: T, TIN, TOUT, start, end
  REAL :: full_sumtime, comp_sumtime, compact_avg, full_avg, setup_time

  ! COMPACTION
  INTEGER, ALLOCATABLE :: REMOVE(:)    ! Vector of species indexes to be removed from full mechanism
  INTEGER :: NRMV                      ! Number of species removed (size(REMOVE))
  
  INTEGER :: idx, S
  INTEGER :: NAVG ! Number of iterations in the timing averaging loop
  INTEGER :: NITR ! Number of integration iterations

  INTEGER, ALLOCATABLE :: rLU_IROW(:), rLU_ICOL(:) ! temporary, for display purposes
  INTEGER              :: iSPC_MAP(NVAR), nonzerocount
  LOGICAL, ALLOCATABLE :: tDO_FUN(:), tDO_SLV(:), tDO_JVS(:)

  REAL(dp)             :: dcdt(NVAR), RxR(NREACT), cinit(NSPEC), lim

  ! Formatting vars
  character(len=20) :: lunz, nv, clunz, cnv

  NITR = 1
  NAVG = 2

  lim = 1e-2

  ALLOCATE(tDO_SLV(NVAR+1))
  ALLOCATE(tDO_FUN(NVAR))
  ALLOCATE(tDO_JVS(LU_NONZERO))

  ALLOCATE(DO_SLV(NVAR+1))  ! Yes/no (1/0), Used to control KppSolve()
  ALLOCATE(DO_FUN(NVAR)  )  ! Yes/no (1/0), Used to control Fun()
  ALLOCATE(DO_JVS(LU_NONZERO)) ! Yes/No (1/0), compute term, used to control Jac_SP()

  call cpu_time(start)

  ! Temporoary settings.
  DO_SLV   = .true. ! Compute all in KppSolve()
  DO_FUN   = .true. ! Compute all in Fun()
  DO_JVS   = .true. ! Compute all in JacSP()
  cNONZERO = 0 ! Initialize number of nonzero elements in reduced mechanism

  tDO_SLV  = .true.
  tDO_FUN  = .true.
  tDO_JVS  = .true.

  iSPC_MAP = 0

  RxR   = 0._dp
  Cinit = 0._dp

  call set_quants_uppertrop(dcdt,Cinit,RxR)

  !where(R.gt.1e-7) R=0.

  ! -------------------------------------------------------------------------- !
  ! 1. Reconstruct the sparse data for a reduced mechanism
  ! e.g. compact the Jacobian

  NRMV     = 0 ! Remove 1 species for testing. This would be determined online.
  do i=1,NVAR
     if (abs(dcdt(i)).le.lim) NRMV=NRMV+1
!     if (abs(dcdt(i)).eq.0.) NRMV=NRMV+1
  ENDDO
  rNVAR    = NVAR-NRMV ! Number of active species in the reduced mechanism

  ! ALLOCATE
  ALLOCATE(REMOVE(NRMV))

  idx = 1
  write(*,*) 'THESE SPECIES ARE INACTIVE: '
  do i=1,NVAR
     if (abs(dcdt(i)).le.lim) then
        write(*,*) 'SPC: ', SPC_NAMES(i), dcdt(i)
        REMOVE(idx)=i
        idx=idx+1
     endif
  ENDDO
  write(*,*) 'NVAR: ', NVAR, '  rNVAR: ', rNVAR, '  ', dble(rNVAR)/dble(NVAR)


  ! -- remove row & column
  ! -- -- Which species are zeroed?


  ! -- DO_SLV, DO_FUN, and DO_JVS will not change size (remain NVAR & NONZERO)
  !    But the appropriate elements are set to zero, so the appropriate terms
  !    in KppSolve(), Fun() and Jac_SP() are not computed

  tDO_SLV(REMOVE(:)) = .false.
  tDO_FUN(REMOVE(:)) = .false.
!  DO_SLV(REMOVE(:)) = .false.
!  DO_FUN(REMOVE(:)) = .false.

  ! -- -- Loop through vectors
  ! -- -- -- Count the number of nonzero elements in the reduced
  ! Jacobian
  nonzerocount = 0
  DO i = 1,LU_NONZERO !NONZERO is the size of LU_ICOL & LU_IROW
!  DO S = 1,NRMV
     IF (.not.(ANY(REMOVE.eq.LU_IROW(i)).or.ANY(REMOVE.eq.LU_ICOL(i)))) THEN
        nonzerocount = nonzerocount+1
     ENDIF
!  ENDDO
  ENDDO

  cNONZERO = nonzerocount

  ! -- -- -- Allocate the new Jacobian elements based on new non-zero elements
  ALLOCATE(cLU_IROW(cNONZERO))
  ALLOCATE(cLU_ICOL(cNONZERO))
  ALLOCATE(JVS_MAP(cNONZERO))
  ALLOCATE(cLU_CROW(rNVAR+1))
  ALLOCATE(cLU_DIAG(rNVAR+1))
  ALLOCATE(SPC_MAP(rNVAR))

  ALLOCATE(rLU_IROW(LU_NONZERO))
  ALLOCATE(rLU_ICOL(LU_NONZERO))
  rLU_IROW = -999
  rLU_ICOL = -999

  ! -- -- Set up SPC_MAP() to map full-mech Fcn to compressed Fcn
  S = 1 ! rNVAR index
  DO i=1,NVAR
     IF (tDO_FUN(i)) THEN
        SPC_MAP(S) = i
        iSPC_MAP(i) = S
        S=S+1
     ENDIF
  ENDDO

  ! -- -- -- recompute LU_IROW
  ! -- -- -- recompute LU_ICOL
  idx = 0
  DO i = 1,LU_NONZERO !NONZERO is the size of LU_ICOL & LU_IROW
     ! Remove row & column elements of deactivated species
     ! Populate cLU_IROW & cLU_ICOL
     IF (ANY(REMOVE.eq.LU_IROW(i)).or.ANY(REMOVE.eq.LU_ICOL(i))) THEN
        rLU_ICOL(i) = 0
        rLU_IROW(i) = 0
        ! Toggle DO_JVS()
        tDO_JVS(i) = .false. ! Turn off this term in Jac_SP()
     ELSE
        idx=idx+1
        cLU_IROW(idx) = iSPC_MAP(LU_IROW(i))
        rLU_IROW(i)   = LU_IROW(i)
        cLU_ICOL(idx) = iSPC_MAP(LU_ICOL(i))
        rLU_ICOL(i)   = LU_ICOL(i)
        JVS_MAP(idx)  = i
     ENDIF
  ENDDO

  cLU_CROW(1) = 1 ! 1st index = 1
  cLU_DIAG(1) = 1 ! 1st index = 1

  ! -- -- -- compute cLU_CROW
  S = 2
  DO i = 2,cNONZERO
     IF ((cLU_IROW(i).ne.cLU_IROW(i-1).and.S.le.rNVAR).or.(i.eq.cNONZERO.and.S.le.rNVAR)) THEN
        cLU_CROW(S) = i
        S=S+1
     ENDIF
  ENDDO

  ! -- -- -- compute cLU_DIAG
  S = 1
  DO i = 2,cNONZERO
     IF (cLU_IROW(i).eq.cLU_ICOL(i)) THEN
        S=S+1
        cLU_DIAG(S) = i 
     ENDIF
  ENDDO

  cLU_CROW(rNVAR+1) = cNONZERO+1
  cLU_DIAG(rNVAR+1) = cLU_DIAG(rNVAR)+1

  call cpu_time(end)
  setup_time = end-start

  ! Initialize formatting values
  write(lunz,*) LU_NONZERO
  write(nv,*) NVAR+1
  write(clunz,*) cNONZERO
  write(cnv,*) rNVAR+1

!  write(*,*) 'cNONZERO: ', cNONZERO
  write(*,*) ' '
  write(*,*) 'Species Info:'
  write(*,*) '---------------'
  do i=1,nvar
!     write(*,'(a,i3)') " Species "//trim(spc_names(i))//" has index: ", i
  enddo
!>>  write(*,*) '---------------'
!>>  write(*,*) ' '
!>>  write(*,*) 'Full-mech sparse data: '
!>>  write(*,'(a,'//lunz//'i4)') ' LU_IROW:  ', LU_IROW
!>>  write(*,'(a,'//lunz//'i4)') ' LU_ICOL:  ', LU_ICOL
!>>  write(*,'(a,'//nv//'i4)') ' LU_CROW:  ', LU_CROW
!>>  write(*,'(a,'//nv//'i4)') ' LU_DIAG:  ', LU_DIAG
!>>  write(*,*) '---------------'
!>>  write(*,*) ' '
!>>  write(*,*) 'Reduced-mech, uncompacted sparse data: '
!  write(*,*) '-- removes species ' // SPC_NAMES(REMOVE(:))! // ' with index ', REMOVE(:)
!>>  write(*,'(a,'//lunz//'i4)') ' rLU_IROW: ', rLU_IROW
!>>  write(*,'(a,'//lunz//'i4)') ' rLU_ICOL: ', rLU_ICOL
!>>  write(*,*) '---------------'
!>>  write(*,*) ' '
!>>  write(*,*) 'Compacted sparse data: '
!>>  write(*,'(a,'//clunz//'i4)') ' cLU_IROW: ', cLU_IROW
!>>  write(*,'(a,'//clunz//'i4)') ' cLU_ICOL: ', cLU_ICOL
!>>  write(*,'(a,'//cnv//'i4)') ' cLU_CROW: ', cLU_CROW
!>>  write(*,'(a,'//cnv//'i4)') ' cLU_DIAG: ', cLU_DIAG
!>>  write(*,*) '---------------'
!>>  write(*,*) ' '
!>>!  write(*,*) 'JVS_MAP ensures that the right JVS values are indexed in the integration'
!>>  write(*,'(a,'//clunz//'i4)') ' JVS_MAP:  ', JVS_MAP
!>>!  write(*,*) ' '
!>>!  write(*,*) 'SPC_MAP ensures that the right species values are indexed in the integration'
!>>  write(*,'(a,'//nv//'i4)') ' SPC_MAP:  ', SPC_MAP
!>>!  write(*,*) ' '
!>>!  write(*,*) 'DO_SLV controls the terms that will be computed in KppSolve(): 1=compute, 0=skip'
!>>  write(*,'(a,'//nv//'l4)') ' DO_SLV:   ', tDO_SLV
!>>!  write(*,*) ' '
!>>!  write(*,*) 'DO_FUN controls the terms that will be computed in Fun(): 1=compute, 0=skip'
!>>  write(*,'(a,'//nv//'l4)') ' DO_FUN:   ', tDO_FUN
!>>!  write(*,*) ' '
!>>!  write(*,*) 'DO_JVS controls the terms that will be computed in Jac_SP(): 1=compute, 0=skip'
!>>  write(*,'(a,'//lunz//'l4)') ' DO_JVS:   ', tDO_JVS
  write(*,*) ' '
  write(*,*) '  Setup time:', setup_time

  ! -------------------------------------------------------------------------- !
  ! 2. Run the full mechanism

  ! Make sure everything is calculated. This should be automatic if not running
  ! a compacted mech. (This is currently also done above, but repeated here for
  ! safety. -- MSL
  DO_FUN = .true.
  DO_SLV = .true.
  DO_JVS = .true.

  call fullmech(NVAR,Cinit,RxR)

  call fullmech(NVAR,Cinit,RxR)
!  call compactedmech()
!  DO i=1,NVAR
!     write(*,*) SPC_NAMES(i), C(i), Cinit(i)
     write(*,*) SPC_NAMES(ind_O3), C(ind_O3), Cinit(ind_O3)
     write(*,*) SPC_NAMES(ind_OH), C(ind_OH), Cinit(ind_OH)
     write(*,*) SPC_NAMES(ind_NO), C(ind_NO), Cinit(ind_NO)
!  END DO

  ! -------------------------------------------------------------------------- !
  ! 3. Run the compacted mechanism
  ! - Need to somehow pass the compacted vectors to KPP
  ! - Should be pretty straight-forward
  ! - Will need to respond to new 'parameters' cNVAR, cNONZERO
  ! - the compacted mechanism will still process the full species vector, 
  !   C(NVAR), not c(rNVAR), only dC/dt of inactive species is zero

  ! OK, now turn on the comp controls

!  call set_quants_uppertrop(dcdt,Cinit,RxR)
  DO_SLV = tDO_SLV
  DO_FUN = tDO_FUN
  DO_JVS = tDO_JVS

!  call fullmech()
  call compactedmech(NVAR,Cinit,RxR)
!  DO i=1,NVAR
     write(*,*) SPC_NAMES(ind_O3), C(ind_O3), Cinit(ind_O3)
     write(*,*) SPC_NAMES(ind_OH), C(ind_OH), Cinit(ind_OH)
     write(*,*) SPC_NAMES(ind_NO), C(ind_NO), Cinit(ind_NO)
!  END DO

  ! -------------------------------------------------------------------------- !
  ! 4. (optional) Calculate the error norm per Santillana et al. (2010) and
  !    Shen et al. (2020)
  !    --- this is valuable if the species removed are selected per
  !        a reduction criterion. Requires a properly posed mechanism.

  ! -------------------------------------------------------------------------- !

  ! 5. Report timing comparison
  ! XXXXX -- For some reason, the ratio below changes depending on which
  !          integration is called first.

  write(*,'(a,f5.1,a)') ' compact/full: ', 100.*compact_avg/full_avg, "%" 
!  write(*,'(a,f4.1,a)') ' problem size: ', 100.*(rNVAR**2)/(NVAR**2), "%"
!  write(*,'(a,f4.1,a)') ' non-zero elm: ', 100.*(cNONZERO)/(LU_NONZERO), "%"

  DEALLOCATE(tDO_SLV)
  DEALLOCATE(tDO_FUN)
  DEALLOCATE(tDO_JVS)
  DEALLOCATE(DO_SLV)
  DEALLOCATE(DO_FUN)
  DEALLOCATE(DO_JVS)
  DEALLOCATE(REMOVE)
  DEALLOCATE(cLU_IROW)
  DEALLOCATE(cLU_ICOL)
  DEALLOCATE(JVS_MAP)
  DEALLOCATE(cLU_CROW)
  DEALLOCATE(cLU_DIAG)
  DEALLOCATE(SPC_MAP)
  DEALLOCATE(rLU_IROW)
  DEALLOCATE(rLU_ICOL)

CONTAINS

  subroutine fullmech( Nv, c_i, RxR )
    USE GCKPP_INTEGRATOR
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE
    IMPLICIT NONE

    integer,  intent(in) :: Nv
    real(dp), intent(in) :: c_i(Nv), RxR(NREACT)
    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0                 ! Rosenbrock output 
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
    TOUT = T + 600._dp
    TEMP = 298.
    
    R(1:NREACT) = RxR(1:NREACT)

    write(*,*) ' '
    write(*,*) 'Running the full mechanism'
    
    full_avg     = 0.
    full_sumtime = 0.
    start        = 0.
    end          = 0.

    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG
       call cpu_time(start)
       DO N=1,NITR
          C(1:NSPEC) = C_i(1:NSPEC)
          ! Initialize
          call Initialize()
          VAR(1:NVAR) = C(1:NVAR)
          FIX(1:NFIX) = C(NVAR+1:NSPEC)
          ! Set RCONST
          call Update_RCONST()
          ! Integrate
          CALL Integrate( TIN,    TOUT,    ICNTRL,      &
               RCNTRL, ISTATUS, RSTATE, IERR )
          C(1:NVAR)   = VAR(:)
!          write(*,*) ISTATUS
       ENDDO
       call cpu_time(end)
       full_sumtime = full_sumtime+end-start
    ENDDO
    full_avg = full_sumtime/dble(NAVG)
    write(*,*) "Average integration time: ", full_avg
    write(*,*) '---------------'
 
    return
 end subroutine fullmech

  subroutine compactedmech( Nv, c_i, RxR )
    USE compact_Integrator
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE

    IMPLICIT NONE

    integer,  intent(in) :: Nv
    real(dp), intent(in) :: c_i(NSPEC), RxR(NREACT)
    
    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0                 ! Rosenbrock output 
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
    TOUT = T + 600._dp
    TEMP = 298.
    
    R = RxR

    write(*,*) ' '
    write(*,*) 'Running the reduced mechanism'
    
    compact_avg  = 0.
    comp_sumtime = 0.
    start        = 0.
    end          = 0.

    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG ! Iterate to generate average comp time
       call cpu_time(start)
       DO N=1,NITR
          C(1:NSPEC) = C_i(1:NSPEC)
          ! Initialize
          call Initialize()
          VAR = C(1:NVAR)
          FIX = C(NVAR+1:NSPEC)
          ! Set RCONST
          call Update_RCONST()
          ! Integrate
          CALL cIntegrate( TIN,    TOUT,    ICNTRL,      &
               RCNTRL, ISTATUS, RSTATE, IERR )
          C(1:NVAR)       = VAR(:)
!          write(*,*) ISTATUS
       ENDDO
       call cpu_time(end)
       comp_sumtime = comp_sumtime+end-start
    ENDDO
    compact_avg = comp_sumtime/dble(NAVG)
    write(*,*) "Average integration time: ", compact_avg
    write(*,*) '---------------'
    
    return
  end subroutine compactedmech

end program main
