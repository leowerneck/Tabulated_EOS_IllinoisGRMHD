
/***********************************************************************************

Helper functions for the interpolation of EOS tables.
The functions are based on  Siegel, Mösta, Desai and Wu 2018, and Schneider, Roberts, Ott, 2017.

***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helpers.h"

void nuc_eos_m_kt1_dpdrhoe_dpderho(const int *restrict n_in,
				     const double *restrict rho, 
				     double *restrict temp,
				     const double *restrict ye,
				     const double *restrict eps,
				     double *restrict dpdrhoe,
				     double *restrict dpderho,
				     const double *restrict prec,
				     int *restrict keyerr,
				     int *restrict anyerr)
{
  const int n = *n_in;
  int keyerr2;

  *anyerr = 0;

  for(int i=0;i<n;i++) {

    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    keyerr2 = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(keyerr2 != 0) {
         *anyerr = 1;
         printf("Inside nuc_eos_m_kt1_dpdrhoe_dpderho kt1 checkbounds,%g, %g, %g",rho[i], temp[i], ye[i]);
         exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));

    int idx[8];
    double delx,dely,delz;
    get_interp_spots(lr,lt,ye[i],&delx,&dely,&delz,idx);
     {
	const int iv = 6;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
     }
     {
	const int iv = 7;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
     }
  }

  return;
}

void nuc_eos_m_kt1_dpdrhot_dedrhot_dpdeps(const int *restrict n_in,
				     const double *restrict rho, 
				     double *restrict temp,
				     const double *restrict ye,
				     double *restrict eps,
				     double *restrict eps2,
				     double *restrict press,
				     double *restrict dpdrho,
				     double *restrict dpdt,
				     double *restrict dedrho,
				     double *restrict dedt,
				     double *restrict dpdeps,
				     const double *restrict prec,
				     int *restrict keyerr,
				     int *restrict anyerr)
{

  const int n = *n_in;
  int keyerr2;

  for(int i=0;i<n;i++) {
    // check if we are fine
    keyerr2 = checkbounds(rho[i], temp[i], ye[i]);
    if(keyerr2 != 0) {
      *anyerr = 1;
      printf("Inside kt1 nuc_eos_m_kt1_dpdrhot_dedrhot_dpdeps checkbounds:%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
  
    int idx[11];
    double delx,dely,delz;
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));

    get_interp_spots_d(lr,lt,ye[i],&delx,&dely,&delz,idx);

    // get prs and eps derivatives by rho and temp

    {
      int iv = 0;
      nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(press[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));
      iv = 1;
      nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(eps[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));
      iv = 7;

      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdeps[i]),iv);

//      dpdeps[i] = exp(dpdeps[i]);
      press[i] = exp(press[i]);
      eps2[i] = exp(eps[i]);
      eps[i] = exp(eps[i])-energy_shift;
    }

  }

  return;
}
void nuc_eos_m_kt1_dpdrhot_dedrhot(const int *restrict n_in,
				     const double *restrict rho, 
				     double *restrict temp,
				     const double *restrict ye,
				     double *restrict eps,
				     double *restrict eps2,
				     double *restrict press,
				     double *restrict dpdrho,
				     double *restrict dpdt,
				     double *restrict dedrho,
				     double *restrict dedt,
				     const double *restrict prec,
				     int *restrict keyerr,
				     int *restrict anyerr)
{

  const int n = *n_in;
  int keyerr2;

  for(int i=0;i<n;i++) {
    // check if we are fine
    keyerr2 = checkbounds(rho[i], temp[i], ye[i]);
    if(keyerr2 != 0) {
      *anyerr = 1;
      printf("Inside kt1 nuc_eos_m_kt1_dpdrhot_dedrhot checkbounds,%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
  
    int idx[11];
    double delx,dely,delz;
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));

    get_interp_spots_d(lr,lt,ye[i],&delx,&dely,&delz,idx);

    // get prs and eps derivatives by rho and temp

    {
      int iv = 0;
      nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(press[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));
      iv = 1;
      nuc_eos_C_linterp_one_d(idx,delx,dely,delz,iv,&(eps[i]), &(dedrho[i]), &(dedt[i]), &(dpdrho[i]), &(dpdt[i]));

      press[i] = exp(press[i]);
      eps2[i] = exp(eps[i]);
      eps[i] = exp(eps[i])-energy_shift;
    }

  }

  return;
}



void nuc_eos_m_kt1_press_eps_mu(const int *restrict n_in,
           const double *restrict rho, 
           const double *restrict temp,
           const double *restrict ye,
           double *restrict eps,
           double *restrict prs,
           double *restrict mu_hat,
           double *restrict mu_e,
           double *restrict mu_p,
           double *restrict mu_n,
           int *restrict keyerr,
           int *restrict anyerr)
{

  const int n = *n_in;
  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    *keyerr = checkbounds(rho[i], temp[i], ye[i]);

    if(*keyerr != 0) {
      *anyerr = 1;

      printf("Inside nuc_eos_m_kt1_press_eps_mu kt1 checkbounds:,%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.

  if(*anyerr) return;
  
  for(int i=0;i<n;i++) {
    //printf("HERE3!\n");
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);
    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);
  
    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }

    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }

    {
      const int iv = 8;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mu_hat[i]),iv);
    }

    {
      const int iv = 9;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mu_e[i]),iv);
      //fprintf(stderr, "mu_e nuc_eos %e\n", (mu_e[i]));
    }

    {
      const int iv = 10;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mu_p[i]),iv);
    }

    {
      const int iv = 11;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(mu_n[i]),iv);
    }
  } 

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
  }

  return;
}

void nuc_eos_m_kt0_press(const int *restrict n_in,
			 const double *restrict rho, 
			 double *restrict temp,
			 const double *restrict ye,
			 const double *restrict eps,
			 double *restrict prs,
			 const double *restrict prec,
			 int *restrict keyerr,
			 int *restrict anyerr)
{

  const int n = *n_in;

  *anyerr = 0;

  for(int i=0;i<n;i++) {
    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    *keyerr = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(*keyerr != 0) {
      *anyerr = 1;
      fprintf(stderr, "Inside nuc_eos_m_kt0_press kt0 checkbounds,%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }


  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;
  // first must find the temperature
  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;
    const double epstot = eps[i]+energy_shift;
    if(epstot>0.0e0) {
      // this is the standard scenario; eps is larger than zero
      // and we can operate with logarithmic tables
      const double lxeps = log(epstot);
      
     // fprintf(stderr, "lr %e,lt %e, ye[i] %e,lxeps %e\n", lr,lt,ye[i],lxeps);
      nuc_eos_findtemp(lr,lt,ye[i],lxeps,*prec, (&ltout), keyerr);
     // fprintf(stderr, "&ltout %e, keyerr %i \n", &ltout,keyerr);

    } else {
      //fprintf(stderr, "keyerr  %i\n", keyerr );
      // will be overwritten further down, only marks error
      *keyerr = 667;
    } // epstot > 0.0
    
    if(*keyerr != 0) {
      // now try negative temperature treatment
      double eps0, eps1;
      int idx[8];
      double delx,dely,delz;



      get_interp_spots_linT_low_eps(lr,temp1,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps1));

      get_interp_spots_linT_low_eps(lr,temp0,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps0));

      //fprintf(stderr, "epstot %e eps0 %e eps1 %e\n", epstot, eps0,eps1 );
      //fprintf(stderr, "temp1 %e temp0 %e \n", temp1, temp0);
      temp[i] = (epstot-eps0) * (temp1-temp0)/(eps1-eps0) + temp0;
      //fprintf(stderr, "temp[i] %e\n", temp[i]);
      // set error codes
      *anyerr = 1;
      *keyerr = 668;

      // get pressure
      {
	const int iv = 0;
	get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
	nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(prs[i]),iv);
      }

    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      {
	const int iv = 0;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
    }
  } // loop i<n
  
  // now get rid of ln:
  for(int i=0;i<n;i++) {
   // fprintf(stderr, "temp %e\n", temp[i]);
   // fprintf(stderr, "prs[i] %e\n", prs[i]);

    prs[i] = exp(prs[i]);
  }

  return;
}


void nuc_eos_m_kt1_press_eps(const int *restrict n_in,
           const double *restrict rho, 
           const double *restrict temp,
           const double *restrict ye,
           double *restrict eps,
           double *restrict prs,
           int *restrict keyerr,
           int *restrict anyerr)
{

  const int n = *n_in;
  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    *keyerr = checkbounds(rho[i], temp[i], ye[i]);

    if(*keyerr != 0) {
      *anyerr = 1;

      printf("Inside kt1 nuc_eos_m_kt1_press_eps  checkbounds:  ,%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.

  if(*anyerr) return;
  
  for(int i=0;i<n;i++) {
    //printf("HERE3!\n");
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);
    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);
  
    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }

    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
  } 

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
  }

  return;
}



void nuc_eos_m_kt0_press_from_entropy(const int *restrict n_in,
       const double *restrict rho, 
       double *restrict temp,
       const double *restrict ye,
       const double *restrict entr,
       double *restrict prs,
       const double *restrict prec,
       int *restrict keyerr,
       int *restrict anyerr)
{

  const int n = *n_in;
  //double ent_test[n];

  *anyerr = 0;

  for(int i=0;i<n;i++) {
    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    *keyerr = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(*keyerr != 0) {
      *anyerr = 1;
      fprintf(stderr, "Inside kt0 nuc_eos_m_kt0_press_from_entropy checkbounds,%g, %g, %g",rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }


  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;
  // first must find the temperature
  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;

      const double xentr = entr[i];
      
      nuc_eos_findtemp_entropy(lr,lt,ye[i],xentr,*prec,(&ltout),keyerr);
     // fprintf(stderr, "&ltout %e, keyerr %i \n", &ltout,keyerr);

    
    if(*keyerr != 0) {
      // now try negative temperature treatment
      double entr0, entr1;
      int idx[8];
      double delx,dely,delz;



      get_interp_spots_linT_low(lr,temp1,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(entr1),2);

      get_interp_spots_linT_low(lr,temp0,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(entr0),2);

      //fprintf(stderr, "epstot %e eps0 %e eps1 %e\n", epstot, eps0,eps1 );
      //fprintf(stderr, "temp1 %e temp0 %e \n", temp1, temp0);
      temp[i] = (xentr-entr0) * (temp1-temp0)/(entr1-entr0) + temp0;
      //fprintf(stderr, "temp[i] %e\n", temp[i]);
      // set error codes
      *anyerr = 1;
      *keyerr = 668;

      // get pressure
      {
  const int iv = 0;
  get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(prs[i]),iv);
      }

    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      {
  const int iv = 0;
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
  //    {
  //nuc_eos_C_linterp_one(idx,delx,dely,delz,&(ent_test[i]),2);
  //    }
    }
  } // loop i<n
  
  // now get rid of ln:
  for(int i=0;i<n;i++) {
   // fprintf(stderr, "temp %e\n", temp[i]);
    prs[i] = exp(prs[i]);
    //fprintf(stderr, "prs[i] %e\n", prs[i]);
    //fprintf(stderr, "ent test %e\n", ent_test[i]);
  }

  return;
}

void nuc_eos_P_from_Enthalpy(double rho, double *temp, double ye,
		     double *enr, double *enr2, double *press,
		     int keytemp,
		     int *keyerr,double rfeps) 

{

  // keyerr codes:
  // 667 -- no temperature found
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  // keytemp codes:
  // 1 -- coming in with temperature
  // 0 -- coming in with eps, need to find temperature
  // 2 -- coming in with spec. entropy, need to find temperature
  //      (not currently implemented)

  int npoints = 1;
  int anyerr = 0;
  *keyerr = 0;

  // set up local vars
 
  double eps2 = 0.0;   
  double press2 = 0.0;   
  double eps = *enr;
  double lr = log(rho);
  double lt = log(*temp);
  double lenthalpy = eps;


  nuc_eos_m_kt1_press_eps(&npoints,&rho,temp,&ye,&eps2,&press2,keyerr,&anyerr);

  //fprintf(stderr, "After nuc_eos_m_kt1_press_eps!\n");

  double enthalpy = eps2+press2/rho;
  enthalpy = enthalpy + enthalpy*0.01;

  if(keytemp == 0) {
    //printf("here Temp: %g Rho: %g Eps: %g Ye: %g\n",*temp,rho,eps,ye);
    double nlt = 0.0;
    nuc_eos_findtemp_enthalpy(lr,lt,ye,lenthalpy,rfeps,&nlt,keyerr);
    if(*keyerr != 0) {
      printf("3 Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",*keyerr,*temp,rho,eps,ye); 
    }
    lt = nlt;

    *temp = exp(lt);
  } else if(keytemp == 1) {
    
  }

  int idx[8];
  double delx,dely,delz,lpress,leps;

  get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&lpress,0);  
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,1);  
  // assign results

  *enr2 = exp(leps) - energy_shift;
  *press= exp(lpress);
  return;
}


void nuc_eos_m_kt1_press_eps_cs2(const int *restrict n_in,
           const double *restrict rho, 
           const double *restrict temp,
           const double *restrict ye,
           double *restrict eps,
           double *restrict prs,
           double *restrict cs2,
           int *restrict keyerr,
           int *restrict anyerr)
{

 // using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    *keyerr = checkbounds(rho[i], temp[i], ye[i]);
    if(*keyerr != 0) {
      *anyerr = 1;
      fprintf(stderr, "Inside kt1 nuc_eos_m_kt1_press_eps_cs2 checkbounds,%g, %g, %g", rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }
  
  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);

    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);
  
    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }
    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
    {
      const int iv = 4;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
    }
  } 

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
    cs2[i] = MAX(cs2[i],0.0e0);
  }

  return;
}

void nuc_eos_m_kt0_press_cs2(const int *restrict n_in,
                         const double *restrict rho,
                         double *restrict temp,
                         const double *restrict ye,
                         const double *restrict eps,
                         double *restrict prs,
                         double *restrict cs2,
                         const double *restrict prec,
                         int *restrict keyerr,
                         int *restrict anyerr)
{


  const int n = *n_in;

  *anyerr = 0;

  for(int i=0;i<n;i++) {

    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    *keyerr = checkbounds_kt0_noTcheck(rho[i], ye[i]);

    if(*keyerr != 0) {
      *anyerr = 1;
      fprintf(stderr, "Inside kt0 nuc_eos_m_kt0_press_cs2 checkbounds,%g, %g, %g",rho[i], temp[i], ye[i]);
      exit(EXIT_FAILURE);
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;
  // first must find the temperature
  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;
    const double epstot = eps[i]+energy_shift;

    if(epstot>0.0e0) {

      // this is the standard scenario; eps is larger than zero
      // and we can operate with logarithmic tables
      const double lxeps = log(epstot);

      nuc_eos_findtemp(lr,lt,ye[i],lxeps,*prec, (&ltout),keyerr);
    } 
    else {
      // will be overwritten further down, only marks error
      *keyerr = 667;
    } // epstot > 0.0

    if(*keyerr != 0) {
      // now try negative temperature treatment
      double eps0, eps1;
      int idx[8];
      double delx,dely,delz;

      get_interp_spots_linT_low_eps(lr,temp1,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps1));

      get_interp_spots_linT_low_eps(lr,temp0,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps0));


      temp[i] = (epstot-eps0) * (temp1-temp0)/(eps1-eps0) + temp0;

      // set error codes
      *anyerr = 1;
      *keyerr = 668;

      // get pressure
      {
        const int iv = 0;
        get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
        nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(prs[i]),iv);
      }
      // get cs2
      {
        const int iv = 4;
        get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
        nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(cs2[i]),iv);
      }
      //fprintf(stderr, "idx %e delx %e dely %e delz %e \n", idx,delx,dely,delz);
      //fprintf(stderr, "lr %e temp[i] %e ye[i] %e &delx %e &dely %e &delz %e idx %e \n", lr,temp[i],ye[i],&delx,&dely,&delz,idx);

    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      {
        const int iv = 0;
        nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
      {
        const int iv = 4;
        nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
      }
      //{
      //  const int iv = 6;
      //  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
     // }
      //{
      //  const int iv = 7;
      //  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
      //}
    }
  }

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    cs2[i] = MAX(cs2[i],0.0e0);
  }

  return;
}

