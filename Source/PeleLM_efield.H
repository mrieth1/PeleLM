// PeleLM efield header. Is included in PeleLM.H if efield requested.
#ifndef _Efield_H_
#define _Efield_H_

// Additional scalar indices
   static int         nE;
   static int         PhiV;
   static int         have_nE;
   static int         have_PhiV;

// MG Poisson solver options	
   static int       			    ef_PoissonMaxIter;
   static int       			    ef_PoissonVerbose;
   static int       			    ef_PoissonMaxOrder;
	static amrex::Real   		 ef_phiV_tol;

// Data
	amrex::MultiFab    kappaSpec_cc;
	amrex::MultiFab    kappaElec_cc;
	amrex::MultiFab    diffElec_cc;

// Functions	
	static void ef_init();

// Define some multifabs for transport properties 
	void ef_define_data();

// SDC predictor : solve for PhiV with 'old' state and compute ef transport properties 
	void ef_advance_setup(amrex::Real time);

// Compute EF related transport properties: Kp_sp, Kp_e, Diff_e
	void ef_calc_transport(amrex::Real time);

// Poisson solve for PhiV with current charge distribution
	void ef_solve_phiv(amrex::Real time);

// Setup PhiV BC for linear Poisson solve
   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc);

#endif /*_Efield_H_*/
