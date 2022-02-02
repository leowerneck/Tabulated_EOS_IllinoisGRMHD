#include "cctk.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>

// EOS_Omni globals
namespace nuc_eos {
  extern double eos_rhomax , eos_rhomin;
  extern double eos_tempmin, eos_tempmax;
  extern double eos_yemin  , eos_yemax;
}
namespace nuc_eos_private {
  extern int nrho,ntemp,nye;
  extern double drho,dtemp,dye;
}

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
  CCTK_VInfo(CCTK_THORNSTRING,"%s: %7.3lf %s",message,t,units);
}

#define CHECK_EOS_OMNI_ERROR(func_name) \
  if( anyerr != 0 || keyerr != 0 ) fprintf(stderr,"(%s) anyerr = %d, keyerr = %d\n",func_name,anyerr,keyerr)

void EOS_Omni_benchmark(CCTK_ARGUMENTS) {

  const double lrmin = log(nuc_eos::eos_rhomin);
  const double ltmin = log(nuc_eos::eos_tempmin);

  FILE *fp1 = fopen("EOS_Omni_press_known_T.dat","w");
  FILE *fp2 = fopen("EOS_Omni_press_unknown_T.dat","w");
  FILE *fp3 = fopen("EOS_Omni_short_known_T.dat","w");
  FILE *fp4 = fopen("EOS_Omni_short_unknown_T.dat","w");
  FILE *fp5 = fopen("EOS_Omni_full_known_T.dat","w");
  FILE *fp6 = fopen("EOS_Omni_full_unknown_T.dat","w");
  double time1,time2,time3,time4,time5,time6;
  struct timespec timer_aux;
  time1 = time2 = time3 = time4 = time5 = time6 = 0.0;
  for(int j=0;j<nuc_eos_private::ntemp;j++) {
    double Y_e = 0.5;
    double T   = exp(ltmin + nuc_eos_private::dtemp*j);
    for(int i=0;i<nuc_eos_private::nrho;i++) {
      double rho = exp(lrmin + nuc_eos_private::drho*i);
      double T_guess;
      double P,eps,S,munu,cs2,dedt,dpderho,dpdrhoe;
      double xa,xh,xn,xp,abar,zbar,mue,mun,mup,muhat;
      int keyerr,anyerr;

      // **********************
      // *** EOS_Omni_press ***
      // **********************
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_press(4,1,1e-8,1,&rho,&eps,&T,&Y_e,&P,&keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_press_1");
      time1 += compute_elapsed_time(&timer_aux);
      fprintf(fp1,"%.15e %.15e %.15e %.15e\n",rho,T,P,eps);
      eps *= 1.1;
      T_guess = 1.1*T;
      if( T_guess > nuc_eos::eos_tempmax ) T_guess = 0.9*nuc_eos::eos_tempmax;
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_press(4,0,1e-8,1,&rho,&eps,&T_guess,&Y_e,&P,&keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_press_0");
      time2 += compute_elapsed_time(&timer_aux);
      fprintf(fp2,"%.15e %.15e %.15e %.15e\n",rho,T_guess,P,eps);

      // **********************
      // *** EOS_Omni_short ***
      // **********************
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_short(4,1,1e-8,1,&rho,&eps,&T,&Y_e,&P,
                     &S,&cs2,&dedt,&dpderho,&dpdrhoe,&munu,
                     &keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_short_1");
      time3 += compute_elapsed_time(&timer_aux);
      fprintf(fp3,"%.15e %.15e %.15e\n",rho,T,S);
      eps *= 1.1;
      T_guess = 1.1*T;
      if( T_guess > nuc_eos::eos_tempmax ) T_guess = 0.9*nuc_eos::eos_tempmax;
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_short(4,0,1e-8,1,&rho,&eps,&T_guess,&Y_e,&P,
                     &S,&cs2,&dedt,&dpderho,&dpdrhoe,&munu,
                     &keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_short_0");
      time4 += compute_elapsed_time(&timer_aux);
      fprintf(fp4,"%.15e %.15e %.15e %.15e\n",rho,T_guess,P,eps);

      // *********************
      // *** EOS_Omni_full ***
      // *********************
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_full(4,1,1e-8,1,&rho,&eps,&T,&Y_e,&P,
                    &S,&cs2,&dedt,&dpderho,&dpdrhoe,
                    &xa,&xh,&xn,&xp,&abar,&zbar,&mue,&mun,&mup,&muhat,
                    &keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_full_1");
      time5 += compute_elapsed_time(&timer_aux);
      fprintf(fp5,"%.15e %.15e %.15e %.15e %.15e %.15e\n",rho,T,mue,mun,mup,muhat);
      eps *= 1.1;
      T_guess = 1.1*T;
      if( T_guess > nuc_eos::eos_tempmax ) T_guess = 0.9*nuc_eos::eos_tempmax;
      start_timer(&timer_aux);
      anyerr = keyerr = 0;
      EOS_Omni_full(4,0,1e-8,1,&rho,&eps,&T_guess,&Y_e,&P,
                    &S,&cs2,&dedt,&dpderho,&dpdrhoe,
                    &xa,&xh,&xn,&xp,&abar,&zbar,&mue,&mun,&mup,&muhat,
                    &keyerr,&anyerr);
      // CHECK_EOS_OMNI_ERROR("EOS_Omni_full_0");
      time6 += compute_elapsed_time(&timer_aux);
      fprintf(fp6,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",rho,T_guess,eps,xa,xh,xn,xp);
    }
    fprintf(fp1,"\n");fprintf(fp2,"\n");
    fprintf(fp3,"\n");fprintf(fp4,"\n");
    fprintf(fp5,"\n");fprintf(fp6,"\n");
  }
  fclose(fp1);fclose(fp2);
  fclose(fp3);fclose(fp4);
  fclose(fp5);fclose(fp6);

  CCTK_INFO("Test #1: Compute P and eps using EOS_Omni_press with known T");
  print_time("Test completed in",time1);
  print_time("Avg. exec. time  ",time1/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  CCTK_INFO("Test #2: Compute P using EOS_Omni_press with unknown T");
  print_time("Test completed in",time2);
  print_time("Avg. exec. time  ",time2/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  CCTK_INFO("Test #3: Compute entropy using EOS_Omni_short with known T");
  print_time("Test completed in",time3);
  print_time("Avg. exec. time  ",time3/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  CCTK_INFO("Test #4: Compute pressure and entropy using EOS_Omni_short with unknown T");
  print_time("Test completed in",time4);
  print_time("Avg. exec. time  ",time4/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  CCTK_INFO("Test #5: Compute chemical potentials using EOS_Omni_full with known T");
  print_time("Test completed in",time5);
  print_time("Avg. exec. time  ",time5/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  CCTK_INFO("Test #6: Compute mass fractions using EOS_Omni_full with unknown T");
  print_time("Test completed in",time6);
  print_time("Avg. exec. time  ",time6/((double)(nuc_eos_private::nrho*nuc_eos_private::ntemp)));

  exit(1);

}
