MODULE mod_bcs_defs

#include <AMReX_REAL.H>

   USE mod_ChemDriver_defs, ONLY: maxspec
 
   IMPLICIT NONE
   
   LOGICAL :: bcinit

   INTEGER, PARAMETER :: Nzones = 5

   REAL_T, DIMENSION(Nzones) :: u_bc 
   REAL_T, DIMENSION(Nzones) :: v_bc
   REAL_T, DIMENSION(Nzones) :: w_bc
   REAL_T, DIMENSION(Nzones) :: rho_bc
   REAL_T, DIMENSION(Nzones) :: T_bc
   REAL_T, DIMENSION(0:maxspec-1,Nzones) :: Y_bc
   REAL_T, DIMENSION(Nzones) :: h_bc

#ifdef USE_EFIELD	
   REAL_T, DIMENSION(Nzones) :: ne_bc
   REAL_T, DIMENSION(Nzones) :: phiV_bc
#endif	

   INTEGER, PARAMETER :: MAXPNTS = 50
   REAL_T, DIMENSION(0:MAXPNTS) :: time_points, vel_points, cntl_points

   CHARACTER(LEN=50) :: ac_hist_file
   REAL_T :: tau_control, cfix, coft_old, sest, V_in_old, corr, &
             changeMax_control, tbase_control, dV_control, scale_control, &
             zbase_control, h_control, controlVelMax
   INTEGER navg_pnts

   INTEGER pseudo_gravity

END MODULE mod_bcs_defs       
