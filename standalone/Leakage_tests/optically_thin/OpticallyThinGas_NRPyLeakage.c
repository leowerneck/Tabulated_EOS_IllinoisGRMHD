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
  const REAL tau_nue[2]  = {0.0,0.0};
  const REAL tau_anue[2] = {0.0,0.0};
  const REAL tau_nux[2]  = {0.0,0.0};
  REAL R_source,Q_source,kappa_nue[2],kappa_anue[2],kappa_nux[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(which_constants,eos_params,rho_b,Y_e,T,tau_nue,tau_anue,tau_nux,
                                                       &R_source,&Q_source,kappa_nue,kappa_anue,kappa_nux);

  // fprintf(stderr,"(NRPyLeakage) R_source = %.15e\n",R_source/rho_b);
  // fprintf(stderr,"(NRPyLeakage) Q_source = %.15e\n",-Q_source/NRPyLeakage_harm_units_cgs_to_geom_Q/1.6021766339999996e-06);

  rhs_gfs[Y_E] = R_source/rho_b;
  rhs_gfs[EPS] = Q_source/rho_b;

  // printf("Input prims: %.15e %.15e %.15e\n", rho_b, Y_e, T);
  // printf("taus       : %.15e %.15e %.15e %.15e %.15e %.15e\n", tau_nue[0], tau_nue[1], tau_anue[0], tau_anue[1], tau_nux[0], tau_nux[1]);
  // printf("kappas     : %.15e %.15e %.15e %.15e %.15e %.15e\n", kappa_nue[0], kappa_nue[1], kappa_anue[0], kappa_anue[1], kappa_nux[0], kappa_nux[1]);
  // printf("R_source   : %.15e\n", R_source);
  // printf("Q_source   : %.15e\n", Q_source);
  // printf("RHSs       : %.15e %.15e\n", rhs_gfs[Y_E], rhs_gfs[EPS]);
}

static inline
void rk4_step_ode( const int which_constants,
                   const NRPyEOS_params *restrict eos_params,
                   const REAL dt,
                   const REAL rho_b,
                   REAL *restrict gfs,
                   REAL *restrict T ) {

  // RK4 (no explicit time dependence on rhs):
  //
  // k1 = dt rhs(y(t))
  // k2 = dt rhs(y(t)+k1/2)
  // k3 = dt rhs(y(t)+k2/2)
  // k4 = dt rhs(y(t)+k3)
  //
  // y(t+dt) = y(t) + (1/6)( k1 + 2(k2 + k3) + k4 )
  REAL k1[2],k2[2],k3[2],k4[2],Y_e,P,eps;

  // RK4 - substep 1
  *T = eos_params->eos_tempmax;
  Y_e = gfs[Y_E];
  eps = gfs[EPS];
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,*T,k1);
  // printf("k1         : %.15e %.15e\n", k1[Y_E], k1[EPS]);

  // RK4 - substep 2;
  *T = eos_params->eos_tempmax;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  // printf("Y_e        : %.15e %.15e %.15e (%.15e)\n", Y_e, gfs[Y_E], k1[Y_E], dt);
  // printf("eps        : %.15e %.15e %.15e (%.15e)\n", eps, gfs[EPS], k1[EPS], dt);
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,*T,k2);
  // printf("k2         : %.15e %.15e\n", k2[Y_E], k2[EPS]);

  // RK4 - substep 3;
  *T = eos_params->eos_tempmax;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  // printf("Y_e        : %.15e %.15e %.15e (%.15e)\n", Y_e, gfs[Y_E], k2[Y_E], dt);
  // printf("eps        : %.15e %.15e %.15e (%.15e)\n", eps, gfs[EPS], k2[EPS], dt);
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,*T,k3);
  // printf("k3         : %.15e %.15e\n", k3[Y_E], k3[EPS]);

  // RK4 - substep 4;
  *T = eos_params->eos_tempmax;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  // printf("Y_e        : %.15e %.15e %.15e (%.15e)\n", Y_e, gfs[Y_E], k3[Y_E], dt);
  // printf("eps        : %.15e %.15e %.15e (%.15e)\n", eps, gfs[EPS], k3[EPS], dt);
  NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho_b,Y_e,eps,&P,T);
  rhs(which_constants,eos_params,rho_b,Y_e,eps,*T,k4);
  // printf("k4         : %.15e %.15e\n", k4[Y_E], k4[EPS]);

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

  // printf("dt = %.15e (%.15e)\n", dt, NRPyLeakage_units_cgs_to_geom_T);

  REAL P,eps;
  NRPyEOS_P_and_eps_from_rho_Ye_T(eos_params,initial_rho_b,initial_Y_e,initial_T,&P,&eps);

  REAL gfs[2] = {initial_Y_e,eps};

  FILE *fp = fopen("opticallythingas_semi_analytic_nrpy.txt","w");
  fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],initial_T);
  // printf("%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],initial_T);
  for(int n=0;n<n_steps;n++) {
    REAL T;
    rk4_step_ode(which_constants, eos_params, dt, initial_rho_b, gfs, &T);
    t += dt;
    // NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,initial_rho_b,gfs[Y_E],gfs[EPS],&P,&T);
    // printf("%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],T);
    fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],T);
  }
  fclose(fp);
}
