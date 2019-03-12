MODULE mod_chemdriver_defs

#include <AMReX_REAL.H>   

   IMPLICIT NONE

!  Parameters   
   INTEGER, PARAMETER :: maxreac  = 472
   INTEGER, PARAMETER :: maxspec  = 72
   INTEGER, PARAMETER :: maxelts  = 6
   INTEGER, PARAMETER :: maxthrdb = 20
   INTEGER, PARAMETER :: maxtp    = 3
   INTEGER, PARAMETER :: maxsp    = 12
   INTEGER, PARAMETER :: maxspnml = 16

!  Acutal size
   INTEGER :: Nelt         ! Number of elements
   INTEGER :: Nspec        ! Number of species
   INTEGER :: Nreac        ! Number of reactions
   INTEGER :: Nfit         ! Size of fitting polynomials

!  Integrator
   DOUBLE PRECISION, DIMENSION(maxspec) ::  spec_scalY
   DOUBLE PRECISION :: vode_rtol, vode_atol

   DOUBLE PRECISION ::  thickFacCH

!  SDC additions
   INTEGER :: negative_Y_test
   DOUBLE PRECISION :: rhoH_INIT, T_cell
   DOUBLE PRECISION, DIMENSION(maxspec) :: rhoY_INIT
   DOUBLE PRECISION, DIMENSION(maxspec+1) ::  c_0, c_1
!$omp threadprivate(c_0, c_1, rhoH_INIT, T_cell, rhoY_INIT, negative_Y_test)

!  Components
   INTEGER :: NP, NRHO, NWT, NWTI, NZ, iN2, iE_sp, NEQ, nchemdiag

!  Typical values
   REAL_T :: typVal_Density, typVal_Temp, typVal_RhoH, typVal_Trac, typVal_Y(maxspec), &
             typVal_Vel, typVal_YMIN, typVal_YMAX

#include <cdwrk.H>  

#ifdef USE_EFIELD
   DOUBLE PRECISION, PARAMETER :: Na = 6.022d23                   ! Avogadro's number
   DOUBLE PRECISION, PARAMETER :: CperECharge = 1.60217662d-19    ! Coulomb per charge
   DOUBLE PRECISION, PARAMETER :: e0 = 8.854187817d-12            ! Free space permittivity (C/(V.m))
   DOUBLE PRECISION, PARAMETER :: er = 1.d0                       ! Relative permittivity of air
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: zk
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: invmwt
   INTEGER, DIMENSION(:), ALLOCATABLE :: spec_charge
#endif		

END MODULE mod_chemdriver_defs
