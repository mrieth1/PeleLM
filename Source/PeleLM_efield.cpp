void PeleLM::ef_advance_setup(Real time) {

   ef_solve_phiv(time); 

	ef_calc_transport(time);
}	

void PeleLM::ef_calc_transport(Real time) {
	BL_PROFILE("EF::ef_calc_transport()");

	const TimeLevel whichTime = which_time(State_Type, time);

	BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

//	Get ptr to current version of diff	
	MultiFab&  diff            = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
   const int nGrow      = diff.nGrow();
 
//	A few check on the grow size of EF transport properties
   BL_ASSERT(kappaSpec_cc.nGrow() >= nGrow);
   BL_ASSERT(kappaElec_cc.nGrow() >= nGrow);
   BL_ASSERT(diffElec_cc.nGrow() >= nGrow);

//	FillPatchIterator for State_type variables. Is it necessary ? For specs, it's been done before for 
//	other transport properties
   FillPatchIterator rhoY_fpi(*this,diff,nGrow,time,State_Type,first_spec,nspecies);
   FillPatchIterator T_fpi(*this,diff,nGrow,time,State_Type,Temp,1);
   FillPatchIterator PhiV_fpi(*this,diff,nGrow,time,State_Type,PhiV,1);
   MultiFab& rhoYmf=rhoY_fpi.get_mf();
   MultiFab& Tmf=T_fpi.get_mf();
   MultiFab& PhiVmf=PhiV_fpi.get_mf();

// Call Fortran to effectively compute the transport properties
// Some of the Fortran should probably be in ChemDriver, but keep it here for now ...	
#ifdef _OPENMP
#pragma omp parallel
#endif  
  for (MFIter mfi(rhoYmf,true); mfi.isValid();++mfi)
  {
     const FArrayBox& RhoD = diff[mfi];
     const FArrayBox& RhoYfab = rhoYmf[mfi];
     const FArrayBox& Tfab = Tmf[mfi];
     const FArrayBox& PhiVfab = PhiVmf[mfi];
     FArrayBox& Kpspfab = kappaSpec_cc[mfi];
     FArrayBox& Kpefab = kappaElec_cc[mfi];
     FArrayBox& Diffefab = diffElec_cc[mfi];
     const Box& gbox = mfi.growntilebox();

     ef_spec_mobility(gbox.loVect(),gbox.hiVect(),
                      Tfab.dataPtr(),    ARLIM(Tfab.loVect()),    ARLIM(Tfab.hiVect()),
                      RhoYfab.dataPtr(), ARLIM(RhoYfab.loVect()), ARLIM(RhoYfab.hiVect()),
                      RhoD.dataPtr(),    ARLIM(RhoD.loVect()),    ARLIM(RhoD.hiVect()),
                      Kpspfab.dataPtr(), ARLIM(Kpspfab.loVect()), ARLIM(Kpspfab.hiVect()));

     ef_elec_mobility(gbox.loVect(),gbox.hiVect(),
                      Tfab.dataPtr(),    ARLIM(Tfab.loVect()),    ARLIM(Tfab.hiVect()),
                      RhoYfab.dataPtr(), ARLIM(RhoYfab.loVect()), ARLIM(RhoYfab.hiVect()),
                      PhiVfab.dataPtr(), ARLIM(PhiVfab.loVect()), ARLIM(PhiVfab.hiVect()),
                      Kpefab.dataPtr(),  ARLIM(Kpefab.loVect()),  ARLIM(Kpefab.hiVect()));

     ef_elec_diffusivity(gbox.loVect(),gbox.hiVect(),
                         Tfab.dataPtr(),     ARLIM(Tfab.loVect()),     ARLIM(Tfab.hiVect()),
                         RhoYfab.dataPtr(),  ARLIM(RhoYfab.loVect()),  ARLIM(RhoYfab.hiVect()),
                         PhiVfab.dataPtr(),  ARLIM(PhiVfab.loVect()),  ARLIM(PhiVfab.hiVect()),
                         Kpefab.dataPtr(),   ARLIM(Kpefab.loVect()),   ARLIM(Kpefab.hiVect()),
							    Diffefab.dataPtr(), ARLIM(Diffefab.loVect()), ARLIM(Diffefab.hiVect()));
  }
	showMF("pnp",diffElec_cc,"pnp_De_cc",level);
	showMF("pnp",kappaElec_cc,"pnp_Ke_cc",level);
}

void PeleLM::ef_define_data() {
   kappaSpec_cc.define(grids,dmap,nspecies,1);
   kappaElec_cc.define(grids,dmap,1,1);
   diffElec_cc.define(grids,dmap,1,1);

	pnp_U.define(grids,dmap,2,1);
	pnp_res.define(grids,dmap,2,1);
	pnp_dU.define(grids,dmap,2,1);
	pnp_bgchrg.define(grids,dmap,1,1);

	pnp_Ueff = new MultiFab[BL_SPACEDIM];
}

void PeleLM::ef_solve_phiv(Real time) {

	MultiFab&  S = get_new_data(State_Type);
	MultiFab rhs_poisson(grids,dmap,1,nGrowAdvForcing);

// Use FillPatchIterator (FPI) to update the data in the growth cell and copy back into S
	FillPatchIterator PhiVfpi(*this,S,1,time,State_Type,PhiV,1);
	MultiFab& PhiVmf = PhiVfpi.get_mf();
	MultiFab::Copy(S,PhiVmf,0,PhiV,1,1);
							
// Get alias to PhiV in S
	MultiFab PhiV_alias(S,amrex::make_alias,PhiV,1);

// Get RHS		
#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
   {
      const Box& box = mfi.tilebox();
      FArrayBox& rhs = rhs_poisson[mfi];
      const FArrayBox& rhoY = S[mfi];
      const FArrayBox& ne = S[mfi];
      ef_calc_rhs_poisson(box.loVect(), box.hiVect(),
                         rhs.dataPtr(0),           ARLIM(rhs.loVect()),    ARLIM(rhs.hiVect()),
					          rhoY.dataPtr(first_spec), ARLIM(rhoY.loVect()),   ARLIM(rhoY.hiVect()),
					          ne.dataPtr(nE),           ARLIM(ne.loVect()),     ARLIM(ne.hiVect()));
   }
	showMF("pnp",rhs_poisson,"pnp_rhspoisson",level);

// Set-up solver tolerances
   const Real S_tol     = ef_phiV_tol;
   const Real S_tol_abs = rhs_poisson.norm0() * ef_phiV_tol;

// Set-up Poisson solver
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);

   MLPoisson phiV_poisson({geom}, {grids}, {dmap}, info);

   phiV_poisson.setMaxOrder(ef_PoissonMaxOrder);

// Set-up BC's
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
	ef_set_PoissonBC(mlmg_lobc, mlmg_hibc);
   phiV_poisson.setDomainBC(mlmg_lobc, mlmg_hibc);
   phiV_poisson.setLevelBC(0, &PhiV_alias);

// LinearSolver options
	MLMG mlmg(phiV_poisson);
	mlmg.setMaxIter(ef_PoissonMaxIter);
	mlmg.setMaxFmgIter(20);
	mlmg.setVerbose(ef_PoissonVerbose);
//	mlmg.setBottomVerbose(bottom_verbose);


// Actual solve
	mlmg.solve({&PhiV_alias}, {&rhs_poisson}, S_tol, S_tol_abs);

}

void PeleLM::ef_init() {
   amrex::Print() << " Init EFIELD solve options \n" << "\n";

   PeleLM::nE                        = -1;
   PeleLM::PhiV                      = -1;
   PeleLM::have_nE                   = 0;
   PeleLM::have_PhiV                 = 0;
   PeleLM::ef_phiV_tol               = 2.0e-8;
   PeleLM::ef_PoissonMaxIter			 = 1000;
   PeleLM::ef_PoissonVerbose			 = 1;
   PeleLM::ef_PoissonMaxOrder			 = 4;
	PeleLM::ef_max_NK_ite             = 20;

	ParmParse pp("ef");

	pp.query("MG_itermax",ef_PoissonMaxIter);
	pp.query("MG_verbose",ef_PoissonVerbose);
	pp.query("MG_maxorder",ef_PoissonMaxOrder);
	pp.query("MG_PhiV_tol",ef_phiV_tol);
}


void PeleLM::ef_solve_PNP(Real dt, 
							     Real time, 
							     MultiFab& Dn,
								  MultiFab& Dnp1,
								  MultiFab& Dhat) {

	BL_PROFILE("EF::ef_solve_PNP()");

// Get edge-averaged transport properties	
   FluxBoxes diff_e(this, 1, 0);
   FluxBoxes conv_e(this, 1, 0);
   MultiFab** kappaElec_ec = conv_e.get();
   MultiFab** diffElec_ec = diff_e.get();
   ef_get_edge_transport(kappaElec_ec, diffElec_ec); 

// Copy some stuff ? Need: macvel, Diff_e, Kp_e, SDC_force, Godunov_force	
// Need the old version of nE
// Need to store the GC of the old version. PNP solved only for interior points

// Get PNP components  	
	MultiFab&  S = get_new_data(State_Type);
	MultiFab&  S_old = get_old_data(State_Type);
	MultiFab nE_old_alias(S_old,amrex::make_alias,nE,1);

	MultiFab::Copy(pnp_U, S, nE, 0, 2, 1);
// Get X scale matrix data : TypValue of nE and PhiV

// Compute provisional CD
   ef_bg_chrg(dt, Dn, Dnp1, Dhat);

// Pre-Newton stuff	
// Get a MF for just nE and PhiV ?
// Get the initial residuals
   ef_NL_residual( nE_old_alias, kappaElec_ec, diffElec_ec, dt );
   const Real norm_NL_res0 = ef_NL_res_norm();
// Get the residual scaling
// Check for direct convergence
// Get data for globalization algo: in 1D I call Jac ... not great ...
//
   Real norm_NL_res = norm_NL_res0; 

bool exit_newton = false;
int NK_ite = 0;	
do {

	NK_ite += 1;
	amrex::Print() << " Newton it: " << NK_ite << " residual: " << norm_NL_res << "\n";

// Init Newton update	
   pnp_dU.setVal(0.0); 	

// GMRES 
//   ef_GMRES_solve(); 	

// Linesearch
   Real lambda = 1.0; 

// Update Newton solution
   pnp_dU.mult(lambda,0,2); 
   pnp_U.plus(pnp_dU,0,2,1);
//   ef_NL_residual( nE_old_alias, kappaElec_ec, diffElec_ec, dt );
   norm_NL_res = ef_NL_res_norm();

// Test exit conditions  	
	test_exit_newton(NK_ite, norm_NL_res0, norm_NL_res, exit_newton);

} while( ! exit_newton );

// Post newton stuff

}

Real PeleLM::ef_NL_res_norm() {
	Real norm = pnp_res.norm2();
	return norm;
}

void PeleLM::test_exit_newton(const int NK_ite, 
									   const Real norm_res0,
										Real norm_res,
										bool& exit_newton) {

  const Real tol_Newton = pow(2.0e-16,1.0/3.0);
  Real max_res = pnp_res.norm0();
  if ( max_res <= tol_Newton ) {
	  exit_newton = true;
	  amrex::Print() << " Converged Newton ite of PNP solve \n";
  }

  if ( NK_ite > ef_max_NK_ite ) {
	  exit_newton = true;
	  amrex::Print() << " Max Newton iteration reached for PNP solve \n";
  }

}

void PeleLM::ef_NL_residual(MultiFab&   ne_old,
								    MultiFab*   Ke_ec[BL_SPACEDIM],
	  								 MultiFab*   De_ec[BL_SPACEDIM],
								    Real        dt) {

   MultiFab& I_R = get_new_data(RhoYdot_Type);
	MultiFab I_R_e(I_R,amrex::make_alias,20,1); // TODO : define iE_sp in C++
 
// Use amrex operators to get the RHS
// Diffusion of ne
	MultiFab diff_ne_term(grids, dmap, 1, 0);
   compute_ne_diffusion_term(dt, De_ec, diff_ne_term);

// Convection of ne
	MultiFab conv_ne_term(grids, dmap, 1, 0);
   compute_ne_convection_term(conv_ne_term);
 
// Laplacian of PhiV
	MultiFab laplacian_term(grids, dmap, 1, 0);
   compute_phiV_laplacian_term(dt, laplacian_term);

// Build the non-linear residual	
// res(ne(:)) = dt * ( diff(:) + conv(:) + I_R(:) ) - ( ne(:) - ne_old(:) )
// res(phiv(:)) = \Sum z_k * \tilde Y_k / q_e - ne + Lapl_PhiV
#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter mfi(pnp_res,true); mfi.isValid(); ++mfi)
   {
      const Box& box = mfi.tilebox();
      const FArrayBox& Ufab = pnp_U[mfi];
      const FArrayBox& difffab = diff_ne_term[mfi];
      const FArrayBox& convfab = conv_ne_term[mfi];
      const FArrayBox& laplfab = laplacian_term[mfi];
      const FArrayBox& IRefab = I_R_e[mfi];
      const FArrayBox& neoldfab = ne_old[mfi];
      const FArrayBox& bgchargfab = pnp_bgchrg[mfi];
      FArrayBox& Resfab = pnp_res[mfi];
	   Resfab.copy(difffab,box,0,box,0,1);		// Copy diff term into res for ne
//	   Resfab.plus(convfab,box,box,0,0,1);		// Add diff term into res for ne
//	   Resfab.plus(IRefab,box,box,0,0,1);	   // Add forcing term (I_R) into res for ne
//		Resfab.mult(dt,box,0,1);					// times dt
//		Resfab.minus(Ufab,box,box,0,0,1);		// Substract current ne
//		Resfab.plus(neoldfab,box,box,0,0,1);	// Add old ne --> Done with ne residuals
	   Resfab.copy(bgchargfab,box,0,box,1,1);	// Copy bg charge term / q_E into res for phiV
		Resfab.minus(Ufab,box,box,0,1,1);	   // Substract current ne
		Resfab.plus(laplfab,box,box,0,1,1);	   // Add the phiV laplacian --> Done with phiV residuals
   }

	showMF("pnp",pnp_res,"pnp_res",level);

}

void PeleLM::ef_bg_chrg(Real dt,
							   MultiFab& Dn,
								MultiFab& Dnp1,
								MultiFab& Dhat) {

#ifdef _OPENMP
#pragma omp parallel
#endif
   for (MFIter mfi(pnp_bgchrg,true); mfi.isValid(); ++mfi) {
		const Box& box = mfi.tilebox();
	   const FArrayBox& rhoYoldfab = get_old_data(State_Type)[mfi];
	   const FArrayBox& afab = (*aofs)[mfi];
		const FArrayBox& dnfab = Dn[mfi];   
		const FArrayBox& dnp1fab = Dnp1[mfi];
		const FArrayBox& dhatfab = Dhat[mfi];
	   const FArrayBox& rfab = get_new_data(RhoYdot_Type)[mfi];
	   FArrayBox& bgchrgfab = pnp_bgchrg[mfi];
		ef_calc_chargedist_prov(box.loVect(), box.hiVect(),
				                  rhoYoldfab.dataPtr(first_spec), ARLIM(rhoYoldfab.loVect()), ARLIM(rhoYoldfab.hiVect()), 
										afab.dataPtr(first_spec), ARLIM(afab.loVect()),       ARLIM(afab.hiVect()),
										dnfab.dataPtr(0),         ARLIM(dnfab.loVect()),      ARLIM(dnfab.hiVect()),
										dnp1fab.dataPtr(0),       ARLIM(dnp1fab.loVect()),    ARLIM(dnp1fab.hiVect()),
										dhatfab.dataPtr(0),       ARLIM(dhatfab.loVect()),    ARLIM(dhatfab.hiVect()),
										rfab.dataPtr(0),          ARLIM(rfab.loVect()),       ARLIM(rfab.hiVect()),
										bgchrgfab.dataPtr(),		  ARLIM(bgchrgfab.loVect()),  ARLIM(bgchrgfab.hiVect()),
										&dt);
   }
}

void PeleLM::ef_get_edge_transport(MultiFab* Ke_ec[BL_SPACEDIM],
											  MultiFab* De_ec[BL_SPACEDIM]) {

   for (MFIter mfi(diffElec_cc,true); mfi.isValid(); ++mfi) {
	   const Box& box = mfi.tilebox();

      for (int dir = 0; dir < BL_SPACEDIM; dir++) {
         FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(nE).lo(dir));
	 	   FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(nE).hi(dir));
	 	   const Box& ebox = mfi.nodaltilebox(dir);
	 	   center_to_edge_fancy((diffElec_cc)[mfi],(*De_ec[dir])[mfi],
		                         amrex::grow(box,amrex::BASISV(dir)), ebox, 0, 
		                         0, 1, geom.Domain(), bc_lo, bc_hi);
	 	   center_to_edge_fancy((kappaElec_cc)[mfi],(*Ke_ec[dir])[mfi],
		                         amrex::grow(box,amrex::BASISV(dir)), ebox, 0, 
		                         0, 1, geom.Domain(), bc_lo, bc_hi);
      }
   }
}

void PeleLM::compute_ne_diffusion_term(Real dt,
													MultiFab* De_ec[BL_SPACEDIM],
													MultiFab& diff_ne) {

// Get alias to ne in pnp_U
//	MultiFab nE_alias(pnp_U,amrex::make_alias,0,1);
	const Real time  = state[State_Type].curTime();	// current time
	MultiFab&  S = get_new_data(State_Type);
	FillPatchIterator nEfpi(*this,S,1,time,State_Type,nE,1);
	MultiFab& nEmf = nEfpi.get_mf();
	MultiFab::Copy(S,nEmf,0,nE,1,1);

// Set-up Lapl operator
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   MLABecLaplacian ne_LAPL({geom}, {grids}, {dmap}, info);
   ne_LAPL.setMaxOrder(ef_PoissonMaxOrder);
	  
// BC's	
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
	ef_set_neBC(mlmg_lobc, mlmg_hibc);
	ne_LAPL.setDomainBC(mlmg_lobc, mlmg_hibc);
	{
      MultiFab nE_alias(S,amrex::make_alias,nE,1);
	   ne_LAPL.setLevelBC(0, &nE_alias);
	}

// Coeff's	
   ne_LAPL.setScalars(0.0, 1.0); 
	std::array<const MultiFab*,AMREX_SPACEDIM> bcoeffs{D_DECL(De_ec[0],De_ec[1],De_ec[2])};
	ne_LAPL.setBCoeffs(0, bcoeffs);

// Get divergence using apply function	
   MultiFab nE_alias(pnp_U,amrex::make_alias,0,1);
	FluxBoxes fluxb  (this, 1, 0);								// Flux box ...
	MultiFab **flux    =   fluxb.get();							// ... associated flux MultiFab
	MLMG mlmg(ne_LAPL);
	mlmg.apply({&diff_ne},{&nE_alias});

	diff_ne.mult(-1.0);

//	std::array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(flux[0],flux[1],flux[2])};
//	mlmg.getFluxes({fp},{&nEU_alias});

// Rescale fluxes to get right stuff regardless of cartesian or r-Z
//   for (int d = 0; d < BL_SPACEDIM; ++d) 
//      flux[d]->mult(1.0/(geom.CellSize()[d]));   

// Get flux divergence, scaled by -1
//	flux_divergence(diff_ne,0,flux,0,1,-1);

}

void PeleLM::compute_ne_convection_term(MultiFab& conv_ne) {

	conv_ne.setVal(0.0);

// Build u_eff : u_mac + Ke*grad(\phi) at faces	
   for (int d = 0; d < BL_SPACEDIM; ++d) 
   	pnp_Ueff[d].copy(u_mac[d]);

}

void PeleLM::compute_phiV_laplacian_term(Real dt,
													  MultiFab& lapl_phiV) {


// Get alias to ne in pnp_U
	const Real time  = state[State_Type].curTime();	// current time
	MultiFab&  S = get_new_data(State_Type);
	FillPatchIterator PhiVfpi(*this,S,1,time,State_Type,PhiV,1);
	MultiFab& PhiVmf = PhiVfpi.get_mf();
	MultiFab::Copy(S,PhiVmf,0,PhiV,1,1);

// Get a bunch of stuff necesary to get the fluxes
	FluxBoxes fluxb  (this, 1, 0);								// Flux box ...
	MultiFab **flux     = fluxb.get();							// ... associated flux MultiFab

// Set-up Poisson operator
   LPInfo info;
   info.setAgglomeration(1);
   info.setConsolidation(1);
   info.setMetricTerm(false);
   MLPoisson phiV_poisson({geom}, {grids}, {dmap}, info);
   phiV_poisson.setMaxOrder(ef_PoissonMaxOrder);

// BC's	
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
	ef_set_PoissonBC(mlmg_lobc, mlmg_hibc);
   phiV_poisson.setDomainBC(mlmg_lobc, mlmg_hibc);
	{
      MultiFab PhiV_alias(S,amrex::make_alias,PhiV,1);
      phiV_poisson.setLevelBC(0, &PhiV_alias);
	}

// LinearSolver to get divergence
	MLMG mlmg(phiV_poisson);
	MultiFab PhiV_alias(pnp_U,amrex::make_alias,1,1);
	mlmg.apply({&lapl_phiV},{&PhiV_alias});

//	std::array<MultiFab*,AMREX_SPACEDIM> fp{D_DECL(flux[0],flux[1],flux[2])};
//	mlmg.getFluxes({fp},{&PhiV_alias});

// Rescale fluxes to get right stuff regardless of cartesian or r-Z
//   for (int d = 0; d < BL_SPACEDIM; ++d) 
//       flux[d]->mult(1.0/(dt*geom.CellSize()[d]));   

// Get flux divergence, scaled by -1
//	flux_divergence(lapl_phiV,0,flux,0,1,-1);

}

// Setup BC conditions for diffusion operator on nE. Directly copied from the diffusion one ...
void PeleLM::ef_set_neBC(std::array<LinOpBCType,AMREX_SPACEDIM>& diff_lobc,
                              std::array<LinOpBCType,AMREX_SPACEDIM>& diff_hibc) {

    const BCRec& bc = get_desc_lst()[State_Type].getBC(nE);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geometry::isPeriodic(idim))
        {
            diff_lobc[idim] = diff_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                diff_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      || 
                     pbc == REFLECT_EVEN)
            {
                diff_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                diff_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                diff_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                diff_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      || 
                     pbc == REFLECT_EVEN)
            {
                diff_hibc[idim] = LinOpBCType::Neumann;
				} 
            else if (pbc == REFLECT_ODD)
            {
                diff_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                diff_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}

// Setup BC conditions for linear Poisson solve on PhiV. Directly copied from the diffusion one ...
void PeleLM::ef_set_PoissonBC(std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_lobc,
                              std::array<LinOpBCType,AMREX_SPACEDIM>& mlmg_hibc) {

    const BCRec& bc = get_desc_lst()[State_Type].getBC(PhiV);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (Geometry::isPeriodic(idim))
        {
            mlmg_lobc[idim] = mlmg_hibc[idim] = LinOpBCType::Periodic;
        }
        else
        {
            int pbc = bc.lo(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      || 
                     pbc == REFLECT_EVEN)
            {
                mlmg_lobc[idim] = LinOpBCType::Neumann;
            }
            else if (pbc == REFLECT_ODD)
            {
                mlmg_lobc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_lobc[idim] = LinOpBCType::bogus;
            }

            pbc = bc.hi(idim);
            if (pbc == EXT_DIR)
            {
                mlmg_hibc[idim] = LinOpBCType::Dirichlet;
            }
            else if (pbc == FOEXTRAP      ||
                     pbc == HOEXTRAP      || 
                     pbc == REFLECT_EVEN)
            {
                mlmg_hibc[idim] = LinOpBCType::Neumann;
				} 
            else if (pbc == REFLECT_ODD)
            {
                mlmg_hibc[idim] = LinOpBCType::reflect_odd;
            }
            else
            {
                mlmg_hibc[idim] = LinOpBCType::bogus;
            }
        }
    }
}
