#include <stdio.h>
#include <time.h>

#include "Basic_defines.h"
#include "NRPyEOS.h"

static inline
void start_timer(struct timespec *a) {
  clock_gettime(CLOCK_REALTIME,a);
}

static inline
double compute_elapsed_time(const struct timespec *a) {
  struct timespec b;
  clock_gettime(CLOCK_REALTIME,&b);
  return( b.tv_sec-a->tv_sec + ((double)(b.tv_nsec-a->tv_nsec))*1e-9 );
}

static inline
void print_time(const char *message, const double time) {
  char units[10];
  double t = time;
  if( t > 1 ) {
    sprintf(units,"s");
  }
  else {
    t *= 1e3;
    if( t > 1 ) {
      sprintf(units,"ms");
    }
    else {
      t *= 1e3;
      sprintf(units,"µs");
    }
  }
  fprintf(stderr,"(NRPyEOS - benchmark) %s: %7.3lf %s\n",message,t,units);
}

void NRPyEOS_benchmark(const NRPyEOS_params *restrict eos_params) {

  struct timespec timer_aux;

  const double lrmin = log(eos_params->eos_rhomin);
  const double ltmin = log(eos_params->eos_tempmin);

  FILE *fp1 = fopen("NRPyEOS_P_and_eps_from_rho_Ye_T.dat","w");
  FILE *fp2 = fopen("NRPyEOS_P_and_T_from_rho_Ye_eps.dat","w");
  FILE *fp3 = fopen("NRPyEOS_S_from_rho_Ye_T.dat","w");
  FILE *fp4 = fopen("NRPyEOS_P_S_and_T_from_rho_Ye_eps.dat","w");
  FILE *fp5 = fopen("NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T.dat","w");
  FILE *fp6 = fopen("NRPyEOS_Xa_Xh_Xn_Xp_and_T_from_rho_Ye_eps.dat","w");
  double time1,time2,time3,time4,time5,time6;
  time1 = time2 = time3 = time4 = time5 = time6 = 0.0;
  for(int j=0;j<eos_params->ntemp;j++) {
    double Y_e = 0.5;
    double T   = exp(ltmin + eos_params->dtemp*j);
    for(int i=0;i<eos_params->nrho;i++) {
      double rho = exp(lrmin + eos_params->drho*i);
      double P,eps;

      // Test #1: NRPyEOS_P_and_eps_from_rho_Ye_T
      start_timer(&timer_aux);
      NRPyEOS_P_and_eps_from_rho_Ye_T(eos_params,rho,Y_e,T,&P,&eps);
      time1 += compute_elapsed_time(&timer_aux);
      fprintf(fp1,"%.15e %.15e %.15e %.15e\n",rho,T,P,eps);

      // Test #2: NRPyEOS_P_and_T_from_rho_Ye_eps
      eps *= 1.1;
      T   *= 1.1;
      T = MAX(MIN(T,0.99*eos_params->eos_tempmax),1.01*eos_params->eos_tempmin);
      start_timer(&timer_aux);
      NRPyEOS_P_and_T_from_rho_Ye_eps(eos_params,rho,Y_e,eps,&P,&T);
      time2 += compute_elapsed_time(&timer_aux);
      fprintf(fp2,"%.15e %.15e %.15e %.15e\n",rho,T,P,eps);

      // Test #3: NRPyEOS_S_from_rho_Ye_T
      double S;
      start_timer(&timer_aux);
      NRPyEOS_S_from_rho_Ye_T(eos_params,rho,Y_e,T,&S);
      time3 += compute_elapsed_time(&timer_aux);
      fprintf(fp3,"%.15e %.15e %.15e\n",rho,T,S);

      // Test #4: NRPyEOS_P_S_and_T_from_rho_Ye_eps
      eps *= 1.1;
      T   *= 1.1;
      T = MAX(MIN(T,0.99*eos_params->eos_tempmax),1.01*eos_params->eos_tempmin);
      start_timer(&timer_aux);
      NRPyEOS_P_S_and_T_from_rho_Ye_eps(eos_params,rho,Y_e,eps,&P,&S,&T);
      time4 += compute_elapsed_time(&timer_aux);
      fprintf(fp4,"%.15e %.15e %.15e %.15e\n",rho,T,P,eps);

      // Test #5: NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T
      double mu_e,mu_p,mu_n,muhat;
      start_timer(&timer_aux);
      NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T(eos_params,rho,Y_e,T,&mu_e,&mu_p,&mu_n,&muhat);
      time5 += compute_elapsed_time(&timer_aux);
      fprintf(fp5,"%.15e %.15e %.15e %.15e %.15e %.15e\n",rho,T,mu_e,mu_p,mu_n,muhat);

      // Test #6: NRPyEOS_P_S_and_T_from_rho_Ye_eps
      double X_a,X_h,X_n,X_p;
      eps *= 1.1;
      T   *= 1.1;
      T = MAX(MIN(T,0.99*eos_params->eos_tempmax),1.01*eos_params->eos_tempmin);
      start_timer(&timer_aux);
      NRPyEOS_Xa_Xh_Xn_Xp_and_T_from_rho_Ye_eps(eos_params,rho,Y_e,eps,&X_a,&X_h,&X_n,&X_p,&T);
      time6 += compute_elapsed_time(&timer_aux);
      fprintf(fp4,"%.15e %.15e %.15e %.15e %.15e %.15e\n",rho,T,X_a,X_h,X_n,X_p);
    }
    fprintf(fp1,"\n");
    fprintf(fp2,"\n");
    fprintf(fp3,"\n");
    fprintf(fp4,"\n");
    fprintf(fp5,"\n");
    fprintf(fp6,"\n");
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);

  fprintf(stderr,"(NRPyEOS - benchmark) Test #1: Compute P and eps using NRPyEOS_P_and_eps_from_rho_Ye_T\n");
  print_time("Test completed in",time1);
  print_time("Avg. exec. time  ",time1/((double)(eos_params->nrho*eos_params->ntemp)));

  fprintf(stderr,"(NRPyEOS - benchmark) Test #2: Compute P and T using NRPyEOS_P_and_T_from_rho_Ye_eps\n");
  print_time("Test completed in",time2);
  print_time("Avg. exec. time  ",time2/((double)(eos_params->nrho*eos_params->ntemp)));

  fprintf(stderr,"(NRPyEOS - benchmark) Test #3: Compute P and eps using NRPyEOS_S_from_rho_Ye_T\n");
  print_time("Test completed in",time3);
  print_time("Avg. exec. time  ",time3/((double)(eos_params->nrho*eos_params->ntemp)));

  fprintf(stderr,"(NRPyEOS - benchmark) Test #4: Compute P and T using NRPyEOS_P_S_and_T_from_rho_Ye_eps\n");
  print_time("Test completed in",time4);
  print_time("Avg. exec. time  ",time4/((double)(eos_params->nrho*eos_params->ntemp)));

  fprintf(stderr,"(NRPyEOS - benchmark) Test #5: Compute chemical potentials using NRPyEOS_mue_mup_mun_and_muhat_from_rho_Ye_T\n");
  print_time("Test completed in",time5);
  print_time("Avg. exec. time  ",time5/((double)(eos_params->nrho*eos_params->ntemp)));

  fprintf(stderr,"(NRPyEOS - benchmark) Test #6: Compute mass fractions using NRPyEOS_Xa_Xh_Xn_Xp_and_T_from_rho_Ye_eps\n");
  print_time("Test completed in",time6);
  print_time("Avg. exec. time  ",time6/((double)(eos_params->nrho*eos_params->ntemp)));
  
}
