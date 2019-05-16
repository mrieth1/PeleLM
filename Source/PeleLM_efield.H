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
	static amrex::Real   		 ef_lambda_jfnk;
   static int       			    ef_max_GMRES_rst;
	static amrex::Real   		 ef_GMRES_reltol;
   static int       			    ef_GMRES_size;

// PNP solve options	
   static int						 ef_max_NK_ite;

// Data
	amrex::MultiFab             kappaSpec_cc;
	amrex::MultiFab             kappaElec_cc;
	amrex::MultiFab             diffElec_cc;
   amrex::MultiFab**           kappaElec_ec;  
   amrex::MultiFab**           diffElec_ec;  

	amrex::MultiFab             pnp_dU;
	amrex::MultiFab             pnp_bgchrg;
	amrex::MultiFab             pnp_gdnv;
	amrex::MultiFab             pnp_nE_old;
	amrex::MultiFab*            pnp_Ueff;
	amrex::MultiFab             pnp_refGC;

// PC data	
	amrex::MLABecCecLaplacian   pnp_pc_diff;
   amrex::MLABecLaplacian      pnp_pc_Stilda;
	amrex::MLABecLaplacian	    pnp_pc_drift;

	amrex::Real			          pnp_SFne;
	amrex::Real			          pnp_SFphiV;
	amrex::Real			          pnp_SUne;
	amrex::Real			          pnp_SUphiV;

// Functions	
	static void ef_init();

// Define some multifabs for transport properties 
	void ef_define_data();

// SDC predictor : solve for PhiV with 'old' state and compute ef transport properties 
	void ef_advance_setup(const amrex::Real &time);

// Compute EF related transport properties: Kp_sp, Kp_e, Diff_e
	void ef_calc_transport(const amrex::Real &time);

// Poisson solve for PhiV with current charge distribution
	void ef_solve_phiv(const amrex::Real &time);

// Solve PNP function: return updated PhiV and nE
	void ef_solve_PNP(const amrex::Real      &dt,
							const amrex::Real      &time, 
						   const amrex::MultiFab  &Dn,
						   const amrex::MultiFab  &Dnp1,
						   const amrex::MultiFab  &Dhat);

// PNP NL residual function
	void ef_NL_residual(const amrex::Real        &dt,
							  const amrex::MultiFab    &pnp_U,
							        amrex::MultiFab    &pnp_res,
							  bool update_scaling = false);

// Test PNP NL residual function
	void ef_NL_residual_test(const amrex::Real       &dt,
							       const amrex::MultiFab   &pnp_U,
							             amrex::MultiFab   &pnp_res,
							       bool update_scaling = false);


// GMRES solve function
   void ef_GMRES_solve(const amrex::Real      &dt,
							  const amrex::Real      &norm_U,
							  const amrex::MultiFab  &U,
							  const amrex::MultiFab  &res,
							        amrex::MultiFab  &dU);

// J times v function	
	void ef_JtV(const amrex::Real      &dt,
				   const amrex::Real      &norm_pnp_U,
				   const amrex::MultiFab  &pnp_U,
				   const amrex::MultiFab  &pnp_res,
				   const amrex::MultiFab  &v_in,
				         amrex::MultiFab  &JtV);

// Preconditioner set-up function. Called once per newton ite. TODO: maybe less ?
   void ef_setUpPrecond(const amrex::Real         &dt,
							         amrex::MultiFab**   De_ec); 

// Apply the preconditioner to any vector v.
   void ef_applyPrecond(const amrex::MultiFab  &v,
								      amrex::MultiFab  &Pv); 	

// Residual norm
   amrex::Real ef_NL_norm(const amrex::MultiFab   &pnp_vec);
   void ef_NL_norm(const amrex::MultiFab             &pnp_vec,
						       amrex::Vector<amrex::Real>  &norm);

// PNP provisional charge distribution
	void ef_bg_chrg(const amrex::Real      &dt,
						 const amrex::MultiFab  &Dn,
						 const amrex::MultiFab  &Dnp1,
						 const amrex::MultiFab  &Dhat);

// PNP exit Newton test
	void test_exit_newton(const amrex::MultiFab &pnp_U,
								 const int             &NK_ite,
								 const amrex::Real     &norm0,
								 const amrex::Real     &norm,
								 bool                  & exit_newton);

// Get electron edge transport properties
	void ef_get_edge_transport(amrex::MultiFab** ke_ec,
										amrex::MultiFab** De_ec);

// Compute the electron diffusion term
   void compute_ne_diffusion_term(const amrex::Real       &dt,
											 const amrex::MultiFab   &pnp_U,
											       amrex::MultiFab** De_ec,
											       amrex::MultiFab   &diff_ne_term);

// Compute the electron convection term
   void compute_ne_convection_term(const amrex::Real       &dt,
											  const amrex::MultiFab   &pnp_U,	
											        amrex::MultiFab** Ke_ec,
											        amrex::MultiFab** phiV_fluxes,
											        amrex::MultiFab   &conv_ne_term);

// Compute the phiV Laplacian term and the PhiV 'fluxes'
   void compute_phiV_laplacian_term(const amrex::Real       &dt,
												const amrex::MultiFab   &pnp_U,
											         amrex::MultiFab** phiV_fluxes,
												      amrex::MultiFab   &phiV_laplacian_term);

// Setup ne BC for diffusion operator
   void ef_set_neBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &diff_lobc,
                    std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &diff_hibc);

// Setup PhiV BC for linear Poisson solve
   void ef_set_PoissonBC(std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_lobc,
                         std::array<amrex::LinOpBCType,AMREX_SPACEDIM> &mlmg_hibc);

#endif /*_Efield_H_*/
