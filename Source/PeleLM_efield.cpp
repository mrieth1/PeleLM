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

void PeleLM::ef_solve_phiv(Real time) {

	MultiFab&  S = get_new_data(State_Type);
	MultiFab rhs_poisson(grids,dmap,1,nGrowAdvForcing);

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

	ParmParse pp("ef");

	pp.query("MG_itermax",ef_PoissonMaxIter);
	pp.query("MG_verbose",ef_PoissonVerbose);
	pp.query("MG_maxorder",ef_PoissonMaxOrder);
	pp.query("MG_PhiV_tol",ef_phiV_tol);
}

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
