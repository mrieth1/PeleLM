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

	void ef_solve_phiv(amrex::Real time);

	static void ef_init();

   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc);

#endif /*_Efield_H_*/
