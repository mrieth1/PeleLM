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

// PNP solve options	
   static int						 ef_max_NK_ite;

// Data
	amrex::MultiFab    kappaSpec_cc;
	amrex::MultiFab    kappaElec_cc;
	amrex::MultiFab    diffElec_cc;

	amrex::MultiFab    pnp_X;
	amrex::MultiFab    pnp_res;
	amrex::MultiFab    pnp_bgchrg;

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

// Solve PNP function: return updated PhiV and nE
	void ef_solve_PNP(amrex::Real      dt,
							amrex::Real    time, 
						   amrex::MultiFab& Dn,
						   amrex::MultiFab& Dnp1,
						   amrex::MultiFab& Dhat);

// PNP NL residual function
	void ef_NL_residual(amrex::MultiFab&   ne_old,
							  amrex::MultiFab** ke_ec,
							  amrex::MultiFab** De_ec,
							  amrex::Real        dt);

// PNP provisional charge distribution
	void ef_bg_chrg(amrex::Real      dt,
						 amrex::MultiFab& Dn,
						 amrex::MultiFab& Dnp1,
						 amrex::MultiFab& Dhat);

// PNP exit Newton test
	void test_exit_newton(const int NK_ite, bool& exit_newton);

// Get electron edge transport properties
	void ef_get_edge_transport(amrex::MultiFab** ke_ec,
										amrex::MultiFab** De_ec);

// Compute the electron diffusion term
   void compute_ne_diffusion_term(amrex::Real      dt,
											 amrex::MultiFab** De_ec,
											 amrex::MultiFab& diff_ne_term);

// Compute the electron convection term
   void compute_ne_convection_term(amrex::MultiFab& conv_ne_term);

// Compute the phiV Laplacian term
   void compute_phiV_laplacian_term(amrex::Real      dt,
												amrex::MultiFab& phiV_laplacian_term);

// Setup ne BC for diffusion operator
   void ef_set_neBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& diff_lobc,
                    std::array<amrex::LinOpBCType,AMREX_SPACEDIM>& diff_hibc);
#endif /*_Efield_H_*/
