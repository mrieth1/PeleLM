  void ne_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                const int* dom_lo, const int* dom_hi,
                const amrex::Real* dx, const amrex::Real* grd_lo,
                const amrex::Real* time, const int* bc);

  void phiv_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                  const int* dom_lo, const int* dom_hi,
                  const amrex::Real* dx, const amrex::Real* grd_lo,
                  const amrex::Real* time, const int* bc);

// Compute RHS of phiV Poisson eq: charge distribution over permittivity 
  void ef_calc_rhs_poisson(const int* lo,  const int* hi,
                           amrex::Real* rhs,          ARLIM_P(rhslo),   ARLIM_P(rhshi),
                           const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                           const amrex::Real* ne,     ARLIM_P(nelo),    ARLIM_P(nehi));

// Get species mobility Kappa_sp from temperature and species diffusivity using Einstein relation
  void ef_spec_mobility(const int* lo, const int* hi,
                        const amrex::Real* tempdat, ARLIM_P(tlo),  ARLIM_P(thi),
                        const amrex::Real* rhoYdat, ARLIM_P(rylo), ARLIM_P(ryhi),
                        const amrex::Real* Diffdat, ARLIM_P(slo),  ARLIM_P(shi),
                        amrex::Real* kpspdat      , ARLIM_P(kplo), ARLIM_P(kphi));

// Get electron mobility Kappa_e : cst value for now 
  void ef_elec_mobility(const int* lo, const int* hi,
                        const amrex::Real* tempdat, ARLIM_P(tlo),   ARLIM_P(thi),
                        const amrex::Real* rhoYdat, ARLIM_P(rylo),  ARLIM_P(ryhi),
                        const amrex::Real* phiVdat, ARLIM_P(plo),   ARLIM_P(phi),
                        amrex::Real* kpedat,        ARLIM_P(kpelo), ARLIM_P(kpehi));

// Get electron diffusivity D_e fron Einstein relation for now 
  void ef_elec_diffusivity(const int* lo, const int* hi,
                           const amrex::Real* tempdat, ARLIM_P(tlo),   ARLIM_P(thi),
                           const amrex::Real* rhoYdat, ARLIM_P(rylo),  ARLIM_P(ryhi),
                           const amrex::Real* phiVdat, ARLIM_P(plo),   ARLIM_P(phi),
                           const amrex::Real* kpedat,  ARLIM_P(kpelo), ARLIM_P(kpehi),
                           amrex::Real* dedat,         ARLIM_P(delo),  ARLIM_P(dehi));
