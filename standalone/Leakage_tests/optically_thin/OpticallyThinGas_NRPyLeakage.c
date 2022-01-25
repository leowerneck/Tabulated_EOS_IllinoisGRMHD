#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"
#include "harm_neutrinos.h"

#define Y_E 0
#define EPS 1

static inline
void rhs( const int which_constants,
          const NRPyEOS_params *restrict eos_params,
          const REAL rho_b,
          const REAL Y_e,
          const REAL eps,
          const REAL T,
          REAL *restrict rhs_gfs ) {
  REAL R_source,Q_source;
  NRPyLeakage_compute_GRMHD_source_terms(which_constants,eos_params,rho_b,Y_e,T,0,0,0,&R_source,&Q_source);

  // fprintf(stderr,"(NRPyLeakage) R_source = %.15e\n",R_source/rho_b);
  // fprintf(stderr,"(NRPyLeakage) Q_source = %.15e\n",-Q_source/NRPyLeakage_harm_units_cgs_to_geom_Q/1.6021766339999996e-06);

  rhs_gfs[Y_E] = R_source/rho_b;
  rhs_gfs[EPS] = Q_source/rho_b;
}

static inline
void rk4_step_ode( const int which_constants,
                   const NRPyEOS_params *restrict eos_params,
                   const REAL dt,
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
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,&T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,T,k1);

  // RK4 - substep 2;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,&T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,T,k2);

  // RK4 - substep 3;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,&T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,T,k3);

  // RK4 - substep 4;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,&T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,T,k4);

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
}

void OpticallyThinGas_NRPyLeakage( const int which_constants,
                                   const NRPyEOS_params *restrict eos_params,
                                   const REAL initial_rho_b,
                                   const REAL initial_Y_e,
                                   const REAL initial_T ) {
  const REAL t_final = 0.5*NRPyLeakage_units_cgs_to_geom_T;
  const REAL dt      = 0.001*NRPyLeakage_units_cgs_to_geom_T;
  const int n_steps  = (int)(t_final/dt+0.5);
  REAL t = 0.0;

  REAL P,eps;
  NRPyEOS_P_and_eps_from_rho_Ye_T(eos_params,initial_rho_b,initial_Y_e,initial_T,&P,&eps);

  REAL gfs[2] = {initial_Y_e,eps};

  FILE *fp = fopen("opticallythingas_semi_analytic_nrpy.txt","w");
  fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],initial_T);
  for(int n=0;n<n_steps;n++) {
    rk4_step_ode(which_constants,eos_params,dt,initial_rho_b,gfs);
    t+= dt;
    REAL T = 1.0;
    NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,initial_rho_b,gfs[Y_E],gfs[EPS],&P,&T);
    fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],T);
  }
  fclose(fp);
}
