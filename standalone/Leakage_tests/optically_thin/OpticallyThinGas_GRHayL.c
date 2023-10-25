//#include "Basic_defines.h"
#include "GRHayL.h"
#include "Neutrinos.h"
#include "NRPyLeakage.h"

#define Y_E 0
#define EPS 1

#define MS_TO_CODE (2.030254467280836e+02)
#define CODE_TO_MS (4.925490947641268e-03)

static inline
void rhs( const eos_parameters *restrict eos,
          const double rho_b,
          const double Y_e,
          const double eps,
          const double T,
          double *restrict rhs_gfs ) {
  neutrino_optical_depths tau = {{0,0},{0,0},{0,0}};
  neutrino_opacities kappa;
  double R_source, Q_source;
  NRPyLeakage_compute_neutrino_opacities_and_GRMHD_source_terms(eos, rho_b, Y_e, T,
                                                                &tau, &kappa, &R_source, &Q_source);

  // fprintf(stderr,"(NRPyLeakage) R_source = %.15e\n",R_source/rho_b);
  // fprintf(stderr,"(NRPyLeakage) Q_source = %.15e\n",-Q_source/NRPyLeakage_harm_units_cgs_to_geom_Q/1.6021766339999996e-06);

  rhs_gfs[Y_E] = R_source/rho_b;
  rhs_gfs[EPS] = Q_source/rho_b;
}

static inline
void rk4_step_ode( const eos_parameters *restrict eos,
                   const double dt,
                   const double rho_b,
                   double *restrict gfs ) {

  // RK4 (no explicit time dependence on rhs):
  //
  // k1 = dt rhs(y(t))
  // k2 = dt rhs(y(t)+k1/2)
  // k3 = dt rhs(y(t)+k2/2)
  // k4 = dt rhs(y(t)+k3)
  //
  // y(t+dt) = y(t) + (1/6)( k1 + 2(k2 + k3) + k4 )
  double k1[2],k2[2],k3[2],k4[2],Y_e,P,eps;
  double T = 1.0;

  // RK4 - substep 1
  Y_e = gfs[Y_E];
  eps = gfs[EPS];
  eos->tabulated_compute_P_T_from_eps(eos, rho_b, Y_e, eps, &P, &T);
  rhs(eos, rho_b, Y_e, eps, T, k1);

  // RK4 - substep 2;
  Y_e = gfs[Y_E] + 0.5*dt*k1[Y_E];
  eps = gfs[EPS] + 0.5*dt*k1[EPS];
  eos->tabulated_compute_P_T_from_eps(eos, rho_b, Y_e, eps, &P, &T);
  rhs(eos, rho_b, Y_e, eps, T, k2);

  // RK4 - substep 3;
  Y_e = gfs[Y_E] + 0.5*dt*k2[Y_E];
  eps = gfs[EPS] + 0.5*dt*k2[EPS];
  eos->tabulated_compute_P_T_from_eps(eos, rho_b, Y_e, eps, &P, &T);
  rhs(eos, rho_b, Y_e, eps, T, k3);

  // RK4 - substep 4;
  Y_e = gfs[Y_E] + dt*k3[Y_E];
  eps = gfs[EPS] + dt*k3[EPS];
  eos->tabulated_compute_P_T_from_eps(eos, rho_b, Y_e, eps, &P, &T);
  rhs(eos, rho_b, Y_e, eps, T, k4);

  // RK4 - update step
  for(int i=0;i<2;i++)
    gfs[i] += (dt/6.0)*( k1[i] + 2.0*( k2[i] + k3[i] ) + k4[i] );
}

void OpticallyThinGas_GRHayL( const char *tablepath,
                              const double initial_rho_b,
                              const double initial_Y_e,
                              const double initial_T ) {

  const double W_max     = 10.0; //IGM default: 10
  const double rho_b_atm = 1e-12;
  const double rho_b_min = -1;
  const double rho_b_max = -1;
  const double Y_e_atm   = 0.5;
  const double Y_e_min   = -1;
  const double Y_e_max   = -1;
  const double T_atm     = 1e-2;
  const double T_min     = -1;
  const double T_max     = -1;

  eos_parameters eos;
  initialize_tabulated_eos_functions_and_params(tablepath, W_max,
                                                rho_b_atm, rho_b_min, rho_b_max,
                                                Y_e_atm, Y_e_min, Y_e_max,
                                                T_atm, T_min, T_max, &eos);
  eos.root_finding_precision=1e-12;

  const double t_final = 0.5*MS_TO_CODE;
  const double dt      = 0.001*MS_TO_CODE;
  const int n_steps  = (int)(t_final/dt+0.5);
  double t = 0.0;

  double P,eps;
  eos.tabulated_compute_P_eps_from_T(&eos, initial_rho_b, initial_Y_e, initial_T, &P, &eps);
  printf("From %.15e %.15e %.15e got %.15e %.15e\n",
         initial_rho_b, initial_Y_e, initial_T, P, eps);
  getchar();

  double gfs[2] = {initial_Y_e, eps};

  FILE *fp = fopen("opticallythingas_semi_analytic_grhayl.txt","w");
  fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*NRPyLeakage_units_geom_to_cgs_T,gfs[Y_E],gfs[EPS],initial_T);
  for(int n=0;n<n_steps;n++) {
    rk4_step_ode(&eos, dt, initial_rho_b, gfs);
    t += dt;
    double T = eos.T_max;
    eos.tabulated_compute_P_T_from_eps(&eos, initial_rho_b, gfs[Y_E], gfs[EPS], &P, &T);
    fprintf(fp,"%.15e %.15e %.15e %.15e\n",t*CODE_TO_MS, gfs[Y_E], gfs[EPS], T);
  }
  fclose(fp);
  eos.tabulated_free_memory(&eos);
}
