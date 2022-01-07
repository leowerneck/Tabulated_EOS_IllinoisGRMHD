#include "Basic_defines.hh"
#include "Leakage.hh"
#include "WVU_EOS_Tabulated_headers.hh"

constexpr int Y_E = 0;
constexpr int EPS = 1;

extern "C" {

inline
void rhs( const REAL rho_b, const REAL Y_e, const REAL T, REAL *restrict rhs_gfs ) {
  REAL R_source,Q_source;
  Leakage_compute_GRMHD_source_terms(rho_b,Y_e,T,1e100,1e100,1e100,&R_source,&Q_source);
  rhs_gfs[Y_E] = R_source/rho_b;
  rhs_gfs[EPS] = Q_source/rho_b;
}

inline
void rk4_step_ode( const REAL dt,
                   const REAL rho_b,
                   REAL *restrict gfs ) {

  // RK4 (no explicit time dependence on rhs):
  //
  // k1 = dt rhs(y(t))
  // k2 = dt rhs(y(t)+k1/2)
  // k3 = dt rhs(y(t)+k2/2)
  // k4 = dt rhs(y(t)+k3)
  //
  // y(t+dt) = y(t) + (1/6)( k1 + 2(k2 + k3) + k4 )
  REAL k1[2],k2[2],k3[2],k4[2],Y_e,P,eps;
  REAL T = 1.0;

  // RK4 - substep 1
  Y_e = gfs[Y_E];
  eps = gfs[EPS];
  WVU_EOS_P_and_T_from_rho_Ye_eps(rho_b,Y_e,eps,&P,&T);
  rhs(rho_b,Y_e,T,k1);

  // RK4 - substep 2;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  WVU_EOS_P_and_T_from_rho_Ye_eps(rho_b,Y_e,eps,&P,&T);
  rhs(rho_b,Y_e,T,k2);

  // RK4 - substep 3;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  WVU_EOS_P_and_T_from_rho_Ye_eps(rho_b,Y_e,eps,&P,&T);
  rhs(rho_b,Y_e,T,k3);

  // RK4 - substep 4;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  WVU_EOS_P_and_T_from_rho_Ye_eps(rho_b,Y_e,eps,&P,&T);
  rhs(rho_b,Y_e,T,k4);

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
}

void Leakage_optically_thin_regime_semi_analytic( const REAL initial_rho_b,
                                                  const REAL initial_Y_e,
                                                  const REAL initial_T ) {
  const REAL t_final = 1.0;
  const REAL dt      = 0.05;
  const int n_steps  = (int)(t_final/dt+0.5);
  REAL t = 0.0;

  REAL P,eps;
  WVU_EOS_P_and_eps_from_rho_Ye_T(initial_rho_b,initial_Y_e,initial_T,&P,&eps);

  REAL gfs[2] = {initial_Y_e,eps};

  FILE *fp = fopen("semi_analytic_optically_thin.txt","w");
  fprintf(fp,"%.15e %.15e %.15e %.15e\n",t,gfs[Y_E],gfs[EPS],initial_T);
  for(int n=0;n<n_steps;n++) {
    rk4_step_ode(dt,initial_rho_b,gfs);
    t+= dt;
    REAL T = 1.0;
    WVU_EOS_P_and_T_from_rho_Ye_eps(initial_rho_b,gfs[Y_E],gfs[EPS],&P,&T);
    fprintf(fp,"%.15e %.15e %.15e %.15e\n",t,gfs[Y_E],gfs[EPS],T);
  }
  fclose(fp);
}

} // extern "C"
