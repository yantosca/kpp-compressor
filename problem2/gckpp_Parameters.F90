! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.3.1_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Parameters.F90
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Parameters

  USE gckpp_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 12 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 12 
! NFLUX - Number of Reaction Flux species
  INTEGER, PARAMETER :: NFLUX = 1 
! NFAM - Number of Prod/Loss Families
  INTEGER, PARAMETER :: NFAM = 1 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 9 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 1 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 12 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 13 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 52 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 53 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 13 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_CO2 = 1 
  INTEGER, PARAMETER :: ind_HNO3 = 2 
  INTEGER, PARAMETER :: ind_H2O = 3 
  INTEGER, PARAMETER :: ind_CH4 = 4 
  INTEGER, PARAMETER :: ind_CO = 5 
  INTEGER, PARAMETER :: ind_H2O2 = 6 
  INTEGER, PARAMETER :: ind_O2 = 7 
  INTEGER, PARAMETER :: ind_OH = 8 
  INTEGER, PARAMETER :: ind_HO2 = 9 
  INTEGER, PARAMETER :: ind_O3 = 10 
  INTEGER, PARAMETER :: ind_NO = 11 
  INTEGER, PARAMETER :: ind_NO2 = 12 

! Index declaration for fixed species in C
!   C(ind_spc)


! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)


END MODULE gckpp_Parameters

