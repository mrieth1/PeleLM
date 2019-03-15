#include <unistd.h>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_Extrapolater.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ErrorList.H>
#include <PeleLM.H>
#include <PeleLM_F.H>
#include <Prob_F.H>
#include <ChemDriver_F.H>
#include <DIFFUSION_F.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiGrid.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_SPACE.H>
#include <AMReX_Interpolater.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_Utility.H>

#if defined(BL_USE_NEWMECH) || defined(BL_USE_VELOCITY)
#include <AMReX_DataServices.H>
#include <AMReX_AmrData.H>
#endif

#include <PROB_F.H>
#include <NAVIERSTOKES_F.H>
#include <DERIVE_F.H>

#include <AMReX_buildInfo.H>

using namespace amrex;

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
}

void PeleLM::ef_define_data() {
   kappaSpec_cc.define(grids,dmap,nspecies,1);
   kappaElec_cc.define(grids,dmap,1,1);
   diffElec_cc.define(grids,dmap,1,1);

	pnp_X.define(grids,dmap,2,1);
	pnp_res.define(grids,dmap,2,1);
	pnp_bgchrg.define(grids,dmap,1,1);
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
                         rhs.dataPtr(),            ARLIM(rhs.loVect()),    ARLIM(rhs.hiVect()),
					          rhoY.dataPtr(first_spec), ARLIM(rhoY.loVect()),   ARLIM(rhoY.hiVect()),
					          ne.dataPtr(nE),           ARLIM(ne.loVect()),     ARLIM(ne.hiVect()));
   }

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
   PeleLM::ef_phiV_tol               = 1.0e-12;
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

void PeleLM::ef_solve_PNP(Real dt, 
							     MultiFab& Dn,
								  MultiFab& Dnp1,
								  MultiFab& Dhat) {

	BL_PROFILE("EF::ef_solve_PNP()");

// Get edge-averaged transport properties	
   FluxBoxes diff_e(this, 1, 0);
   FluxBoxes conv_e(this, 1, 0);
   Multifab** kappaElec_ec = conv_e.get();
   Multifab** diffElec_ec = diff_e.get();
   ef_get_edge_transport(diffElec_ec, kappaElec_ec); 

// Copy some stuff ? Need: macvel, Diff_e, Kp_e, SDC_force, Godunov_force	
// Need the old version of nE
// Need to store the GC of the old version. PNP solved only for interior points

// Get PNP components  	
	MultiFab&  S = get_new_data(State_Type);
	MultiFab&  S_old = get_old_data(State_Type);
	MultiFab nE_old_alias(S_old,amrex::make_alias,nE,1);

	MultiFab::Copy(pnp_X, S, nE, 0, 2, 1);
// Get X scale matrix data : TypValue of nE and PhiV

// Compute provisional CD
   ef_bg_chrg(dt, Dn, Dnp1, Dhat);

// Pre-Newton stuff	
// Get a MF for just nE and PhiV ?
// Get the initial residuals
   ef_NL_residual(dt);
// Get the residual scaling
// Check for direct convergence
// Get data for globalization algo: in 1D I call Jac ... not great ...

bool exit_newton = false;
int NK_ite = 0;	
do {

	NK_ite += 1;
	test_exit_newton(NK_ite, exit_newton);

} while( ! exit_newton );

// Post newton stuff

}

void PeleLM::test_exit_newton(const int NK_ite, bool& exit_newton) {

  if ( NK_ite > ef_max_NK_ite ) {
	  exit_newton = true;
	  amrex::Print() << " Max Newton iteration reached for PNP solve \n";
  }

}

void PeleLM::ef_NL_residual(Real dt) {

// Use amrex operators to get the RHS
 
// Diffusion of ne

// Convection of ne
  
// Laplacian of PhiV

//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//   for (MFIter mfi(pnp_X,true); mfi.isValid(); ++mfi)
//   {
//      const Box& box = mfi.tilebox();
//      const FArrayBox& Xfab = pnp_X[mfi];
//      const FArrayBox& Kpefab = kappaElec_cc[mfi];
//      const FArrayBox& Defab = diffElec_cc[mfi];
//      const FArrayBox& bgchrgfab = pnp_bgchrg[mfi];
//      FArrayBox& Resfab = pnp_res[mfi];
//      ef_calc_NL_residual(box.loVect(), box.hiVect(),
//                          Xfab.dataPtr(0),      ARLIM(Xfab.loVect()),       ARLIM(Xfab.hiVect()),
//                          Kpefab.dataPtr(),     ARLIM(Kpefab.loVect()),     ARLIM(Kpefab.hiVect()),
//                          Defab.dataPtr(),      ARLIM(Defab.loVect()),      ARLIM(Defab.hiVect()),
//                          bgchrgfab.dataPtr(),  ARLIM(bgchrgfab.loVect()),  ARLIM(bgchrgfab.hiVect()),
//				              Resfab.dataPtr(0),    ARLIM(Resfab.loVect()),     ARLIM(Resfab.hiVect()),
//								  dt);
//   }
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
				                  rhoYoldfab.dataPtr(0),    ARLIM(rhoYoldfab.loVect()), ARLIM(rhoYoldfab.hiVect()), 
										afab.dataPtr(first_spec), ARLIM(afab.loVect()),       ARLIM(afab.hiVect()),
										dnfab.dataPtr(0),         ARLIM(dnfab.loVect()),      ARLIM(dnfab.hiVect()),
										dnp1fab.dataPtr(0),       ARLIM(dnp1fab.loVect()),    ARLIM(dnp1fab.hiVect()),
										dhatfab.dataPtr(0),       ARLIM(dhatfab.loVect()),    ARLIM(dhatfab.hiVect()),
										rfab.dataPtr(0),          ARLIM(rfab.loVect()),       ARLIM(rfab.hiVect()),
										bgchrgfab.dataPtr(),		  ARLIM(bgchrgfab.loVect()),  ARLIM(bgchrgfab.hiVect()),
										dt);
   }
}

void PeleLM::ef_get_edge_transport(MultiFab* Ke_ec[BL_SPACEDIM],
											  MultiFab* De_ec[BL_SPACEDIM]) {

   for (MFIter mfi(diffElec_cc,true); mfi.isValid(); ++mfi) {
	   const Box& box = mfi.tilebox();

      for (int dir = 0; dir < BL_SPACEDIM; dir++) {
         FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).lo(dir));
	 	   FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).hi(dir));
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
