  void ne_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                const int* dom_lo, const int* dom_hi,
                const amrex::Real* dx, const amrex::Real* grd_lo,
                const amrex::Real* time, const int* bc);

  void phiv_fill (amrex::Real* data, ARLIM_P(data_lo), ARLIM_P(data_hi),
                  const int* dom_lo, const int* dom_hi,
                  const amrex::Real* dx, const amrex::Real* grd_lo,
                  const amrex::Real* time, const int* bc);

  void ef_calc_rhs_poisson(const int* lo,  const int* hi,
                           amrex::Real* rhs,          ARLIM_P(rhslo),   ARLIM_P(rhshi),
                           const amrex::Real* rhoY,   ARLIM_P(rYlo),    ARLIM_P(rYhi),
                           const amrex::Real* ne,     ARLIM_P(nelo),    ARLIM_P(nehi));
