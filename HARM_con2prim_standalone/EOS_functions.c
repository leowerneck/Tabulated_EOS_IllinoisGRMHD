
/***********************************************************************************

This file has the interpolation routines calls for the EOS tables. The functions are based on 
Siegel, Mösta, Desai and Wu 2018. The interpolation driver comes from stellarcollapse 
(Schneider, Roberts, Ott, 2017).


  -- v1.0 :  Written by Ariadna Murguia-Berthier (April 2019) 


**************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.h"



void EOS_press(double rho, double * eps, 
        double * temp, double ye, double * prs);

void EOS_Omni_press(double prec, double rho, double * eps, 
    double * temp, double ye,  double * press);

void EOS_press_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs);

void EOS_press_mu_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs, 
        double * mu_hat,
        double * mu_e, double * mu_p, double * mu_n);

void EOS_Omni_press_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press);

void EOS_Omni_press_mu_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press, 
    double * mu_hat,
    double * mu_e, double * mu_p, double * mu_n);

void EOS_press_cs2(double rho, double * eps, 
        double * temp, double ye, double * prs, double * cs2);

void EOS_Omni_press_cs2(double prec, double rho, double * eps, 
    double * temp, double ye,  double * press, double * cs2);

void EOS_press_cs2_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs, double * cs2);

void EOS_Omni_press_cs2_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press, double * cs2);


void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, 
          double * temp_guess, double ye, double * eps, double * press, double * dPdrho, 
          double * dPdeps, double * entr,const double prec);

void EOS_Omni_dpdrho_dpdeps_hinv(const double rho, const double enth, double * temp, 
       const double ye, double * eps, double * press, double * dpdrho, double * dpdeps, double * entr, 
      int * keyerr, int stepsize);


extern void nuc_eos_m_kt0_press(const int *restrict n_in,
                         const double *restrict rho,
                         double *restrict temp,
                         const double *restrict ye,
                         const double *restrict eps,
                         double *restrict prs,
                         const double *restrict prec,
                         int *restrict keyerr,
                         int *restrict anyerr);



extern void nuc_eos_m_kt0_press_from_entropy(const int *restrict n_in,
       const double *restrict rho, 
       double *restrict temp,
       const double *restrict ye,
       const double *restrict entr,
       double *restrict prs,
       const double *restrict prec,
       int *restrict keyerr,
       int *restrict anyerr);

extern void nuc_eos_m_kt1_press_eps_mu(const int *restrict n_in,
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
           int *restrict anyerr);

extern void nuc_eos_m_kt1_press_eps_cs2(const int *restrict n_in,
           const double *restrict rho, 
           const double *restrict temp,
           const double *restrict ye,
           double *restrict eps,
           double *restrict prs,
           double *restrict cs2,
           int *restrict keyerr,
           int *restrict anyerr);

extern void nuc_eos_m_kt1_dpdrhot_dedrhot_dpdeps(const int *restrict n_in,
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
             int *restrict anyerr);

extern void nuc_eos_m_kt1_dpdrhot_dedrhot(const int *restrict n_in,
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
             int *restrict anyerr);

extern void nuc_eos_P_from_Enthalpy(double rho, double *temp, double ye,
         double *enr, double *enr2, double *press,
         int keytemp,
         int *keyerr,double rfeps); 

void EOS_Omni_press_entropy(double prec, double rho, double * entr, 
    double * temp, double ye,  double * press);


void EOS_press_entropy(double rho, double * entr, 
        double * temp, double ye, double * prs); 




extern void nuc_eos_m_kt1_dpdrhoe_dpderho(const int *restrict n_in,
                                     const double *restrict rho,
                                     double *restrict temp,
                                     const double *restrict ye,
                                     const double *restrict eps,
                                     double *restrict dpdrhoe,
                                     double *restrict dpderho,
                                     const double *restrict prec,
                                     int *restrict keyerr,
                                     int *restrict anyerr);

void nuc_eos_m_kt0_press_cs2(const int *restrict n_in,
                         const double *restrict rho,
                         double *restrict temp,
                         const double *restrict ye,
                         const double *restrict eps,
                         double *restrict prs,
                         double *restrict cs2,
                         const double *restrict prec,
                         int *restrict keyerr,
                         int *restrict anyerr);

void nuc_eos_m_kt1_press_eps(const int *restrict n_in,
           const double *restrict rho, 
           const double *restrict temp,
           const double *restrict ye,
           double *restrict eps,
           double *restrict prs,
           int *restrict keyerr,
           int *restrict anyerr);


void EOS_Omni_dpdrho_dpdt_dedrho_dedt(double rho, double * temp, 
       double ye, double * eps, double * press, double * dpdrho, double * dpdt, double * dedrho, 
       double * dedt, int * keyerr);


void EOS_EP_dEdr_dEdt_dPdr_dPdt(double * x, const double D, 
          double * Eprim, double * Pprim, double * dEdrho, double * dEdt, double * dPdrho, 
          double * dPdt, double ye);

void EP_dEdgamma_dEdW_dEdT(double * Eprim, double * Pprim, double * dEdgamma, double * dEdW, double * dEdT, 
       double * dpdrho, double * dpdT, double * x, const double D, double ye);



/**********************************************************************/
/******************************************************************

  EOS_press():
  Interpolating EOS table routine. Takes the density and internal energy and calls the driver
  that interpolates the table (), giving the pressure. 

******************************************************************/

void EOS_press(double rho, double * eps, 
        double * temp, double ye, double * prs) {

  const double prec = 1.0e-10;
  double xrho, xeps, xtemp, xye, xprs;

  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;

  
  EOS_Omni_press(prec,xrho,&xeps,&xtemp,xye,&xprs);

  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
}


void EOS_Omni_press(double prec, double rho, double * eps, 
    double * temp, double ye,  double * press){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;

    nuc_eos_m_kt0_press(&npoints,&rho,temp,&ye,eps,press,&prec,&keyerr,&anyerr);

}


/**********************************************************************/
/******************************************************************

  EOS_press_entropy():
  Interpolating EOS table routine. Takes the density and entropy and calls the driver
  that interpolates the table (), giving the pressure. 

******************************************************************/

void EOS_press_entropy(double rho, double * entr, 
        double * temp, double ye, double * prs) {

  const double prec = 1.0e-10;
  double xrho, xentr, xtemp, xye, xprs;

  xrho = rho;
  xtemp = *temp;
  xentr = *entr;
  xye = ye;
  xprs = 0.0;

  
  EOS_Omni_press_entropy(prec,xrho,&xentr,&xtemp,xye,&xprs);

  *entr  = xentr;
  *temp = xtemp;
  *prs  = xprs;
}

void EOS_Omni_press_entropy(double prec, double rho, double * entr, 
    double * temp, double ye,  double * press){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;

    nuc_eos_m_kt0_press_from_entropy(&npoints,&rho,temp,&ye,entr,press,&prec,&keyerr,&anyerr);

}

/**********************************************************************/
/******************************************************************

  EOS_press_with_T():
  Interpolating EOS table routine. Takes the density and temperature calls the driver
  that interpolates the table, giving the pressure. 

******************************************************************/



void EOS_press_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs) {

  //const double prec = 1.0e-12;
  double xrho, xeps, xtemp, xye, xprs;

  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;
  
  EOS_Omni_press_with_T(xrho,&xeps,&xtemp,xye,&xprs);

  
  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
}



void EOS_Omni_press_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;

    nuc_eos_m_kt1_press_eps(&npoints,&rho,temp,&ye,eps,press,&keyerr,&anyerr);

}



/**********************************************************************/
/******************************************************************

  EOS_press_with_T():
  Interpolating EOS table routine. Takes the density and temperature calls the driver
  that interpolates the table, giving the pressure, and chemical potentials. 

******************************************************************/



void EOS_press_mu_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs, 
        double * mu_hat,
        double * mu_e, double * mu_p, double * mu_n) {

  //const double prec = 1.0e-12;
  double xrho, xeps, xtemp, xye, xprs, xmu_hat, xmu_e, xmu_p,xmu_n;

  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;
  xmu_hat = 0.0;
  xmu_e = 0.0;
  xmu_p = 0.0;
  xmu_n = 0.0;

  
  EOS_Omni_press_mu_with_T(xrho,&xeps,&xtemp,xye,&xprs,&xmu_hat,&xmu_e,&xmu_p,&xmu_n);
  
  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
  *mu_hat  = xmu_hat;
  *mu_e  = xmu_e;
  *mu_p  = xmu_p;
  *mu_n  = xmu_n;
}



void EOS_Omni_press_mu_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press, 
    double * mu_hat,
    double * mu_e, double * mu_p, double * mu_n){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;

    nuc_eos_m_kt1_press_eps_mu(&npoints,&rho,temp,&ye,eps,press,mu_hat,mu_e,mu_p,mu_n,&keyerr,&anyerr);

}


/**********************************************************************/
/******************************************************************

  EOS_press_cs2():
  Interpolating EOS table routine. Takes the density and internal calls the driver
  that interpolates the table (), giving the pressure and soundspeed. 

******************************************************************/


void EOS_press_cs2(double rho, double * eps, 
        double * temp, double ye, double * prs, double * cs2) {

  const double prec = 1.0e-12;
  double xrho, xeps, xtemp, xye, xprs, xcs2;


  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;
  xcs2=0.;
  

  EOS_Omni_press_cs2(prec,xrho,&xeps,&xtemp,xye,&xprs, &xcs2);



  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
  *cs2  = xcs2;
}


void EOS_Omni_press_cs2(double prec, double rho, double * eps, 
    double * temp, double ye,  double * press, double * cs2){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;
    

    nuc_eos_m_kt0_press_cs2(&npoints,&rho,temp,&ye,eps,press,cs2,&prec,&keyerr,&anyerr);

    
}

/**********************************************************************/
/******************************************************************

  EOS_press_cs2_with_T():
  Interpolating EOS table routine. Takes the density and temperature calls the driver
  that interpolates the table, giving the pressure, and cs2. 

******************************************************************/



void EOS_press_cs2_with_T(double rho, double * eps, 
        double * temp, double ye, double * prs, double * cs2) {

  //const double prec = 1.0e-12;
  double xrho, xeps, xtemp, xye, xprs,xcs2;

  xrho = rho;
  xtemp = *temp;
  xeps = *eps;
  xye = ye;
  xprs = 0.0;
  xcs2=0.0;
  
  EOS_Omni_press_cs2_with_T(xrho,&xeps,&xtemp,xye,&xprs,&xcs2);

  
  *eps  = xeps;
  *temp = xtemp;
  *prs  = xprs;
  *cs2  = xcs2;
}



void EOS_Omni_press_cs2_with_T(double rho, double * eps, 
    double * temp, double ye,  double * press, double * cs2){

    int npoints = 1;
    int keyerr = 0;
    int anyerr = 0;

    nuc_eos_m_kt1_press_eps_cs2(&npoints,&rho,temp,&ye,eps,press,cs2,&keyerr,&anyerr);

}

/**********************************************************************/
/******************************************************************

  EOS_P_from_hrho_dPdrho_dPdeps():
  Performs inversion from specific enthalpy to temperature and compute partial derivatives 
   of pressure wrt density and specific internal energy

******************************************************************/



void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, 
          double * temp_guess, double ye, double * eps, double * press, double * dPdrho, 
          double * dPdeps, double * entr,const double prec){
  

  const int step_size = 1;
  int keyerr;
  double xtemp,xeps,xprs,xdpdrho,xdpdeps,xentr;
  const double xenth = enth;
  const double xrho = rho;

 

  keyerr = 0;
  xtemp = *temp_guess;
  xeps = 0.0;
  xprs = 0.0;
  xdpdeps = 0.0;
  xdpdrho = 0.0;
  xentr = 0.0;

  

  EOS_Omni_dpdrho_dpdeps_hinv(xrho, xenth, &xtemp, ye, &xeps, &xprs, &xdpdrho, 
        &xdpdeps, &xentr, &keyerr, step_size);
  //fprintf(stderr, "After EOS_Omni_dpdrho_dpdeps_hinv\n");
    
  *temp_guess = xtemp; 
  *dPdeps = xdpdeps;
  *dPdrho = xdpdrho;
  *press  = xprs;
  *eps = xeps;
  *entr = xentr;
  
}


void EOS_Omni_dpdrho_dpdeps_hinv(const double rho, const double enth, double * temp, 
       const double ye, double * eps, double * press, double * dpdrho, double * dpdeps, double * entr, 
      int * keyerr, int stepsize) {



    double entm1;
    const double prec = 1.0e-10;
    const int keytemp = 0; 
    entm1 = enth-1.0;
  
    nuc_eos_P_from_Enthalpy(rho, temp, ye, &entm1, eps, press,
        keytemp, keyerr, prec);

    int anyerr = 0; 
    int npoints =1;  

    nuc_eos_m_kt1_dpdrhoe_dpderho(&npoints,&rho,temp,&ye,eps,
                    dpdrho,dpdeps,&prec,keyerr,&anyerr);

    
    // modify EOS call to also output entropy set to zero for now
    *entr = 0.0;


}



/**********************************************************************/
/******************************************************************

  EP_dEdgamma_dEdW_dEdT():
  Compute partial derivatives of specific internal energy and pressure with respect
  to W and gamma
  Note: E = eps - eps_EOS

******************************************************************/


void EP_dEdgamma_dEdW_dEdT(double * Eprim, double * Pprim, double * dEdgamma, double * dEdW, double * dEdT, 
       double * dpdrho, double * dpdT, double * x, const double D, double ye){

  
  double gamma = x[0];
  double W = x[1];
  double T = x[2];
  
  double epsEOS,pEOS,depsEOSdrho,dpEOSdrho,depsEOSdt,dpEOSdt;
  epsEOS = 0.0;
  pEOS = 0.0;
  depsEOSdrho = 0.0;
  dpEOSdrho=0.0;
  dpEOSdt = 0.0;
  
  /* need partial derivatives of specific internal energy and pressure wrt density and 
   * temperature. Those need to be based on primitives computed from Newton-Raphson state
   * vector x and conservatives
   */

  EOS_EP_dEdr_dEdt_dPdr_dPdt(x,D,&epsEOS,&pEOS,&depsEOSdrho,&depsEOSdt,&dpEOSdrho,&dpEOSdt,ye);

  // Further partial derivatives
  double drhodgamma = -D/(gamma*gamma);
  double depsdgamma = -W/(D*gamma*gamma)-pEOS/D + dpEOSdrho/gamma;
   double depsdP = -gamma/D;
  double drhodW = 0.0;
  double depsdW = 1.0/(D*gamma);

  *Eprim = epsEOS;
  *Pprim = pEOS;
  *dEdgamma = depsdgamma - depsEOSdrho*drhodgamma;
  *dEdW = depsdW;
  *dEdT = depsdP*dpEOSdt - depsEOSdt;
  *dpdrho = dpEOSdrho;
  *dpdT = dpEOSdt;

}


/**********************************************************************/
/******************************************************************

  EP_dEdgamma_dEdW_dEdT():
  Compute partial derivatives of specific internal energy and pressure with respect
  to density and temperature, based on primitives computed from Newton-Raphson state
   vector x and conservatives

******************************************************************/


void EOS_EP_dEdr_dEdt_dPdr_dPdt(double * x, const double D, 
          double * Eprim, double * Pprim, double * dEdrho, double * dEdt, double * dPdrho, 
          double * dPdt, double ye){
  
  
  const double gamma = x[0];
  const double W = x[1];
  const double T = x[2];
 
 
  //const int keytemp = 1;
  const double prec = 1e-10;
  int keyerr;
  double xrho,xtemp,xeps,xye,xprs,xdedrho,xdpdrho,xdedt,xdpdt;

  
  keyerr = 0;
  
  xrho = D/gamma;
  xtemp = T;
  xeps = 0.0;
  xye =ye;
  xprs = 0.0;
  
  xdedrho = 0.0;
  xdpdrho = 0.0;
  xdedt = 0.0;
  xdpdt = 0.0;
  
  EOS_Omni_dpdrho_dpdt_dedrho_dedt(xrho, &xtemp, xye, &xeps, &xprs, &xdpdrho, 
        &xdpdt, &xdedrho, &xdedt, &keyerr);

    
  *dEdrho = xdedrho;
  *dEdt   = xdedt;
  *dPdrho = xdpdrho;
  *dPdt   = xdpdt;
  *Eprim  = xeps;
  *Pprim  = xprs;

  
}


void EOS_Omni_dpdrho_dpdt_dedrho_dedt(double rho, double * temp, 
       double ye, double * eps, double * press, double * dpdrho, double * dpdt, double * dedrho, 
       double * dedt, int * keyerr)
{

// nuc EOS
  
    double xtemp,xprs;
    const double prec = 1.0e-10;
    double xeps2 = 0.0; 
    int anyerr = 0; 

    xtemp = *temp;
    //xeps = *eps * inv_eps_gf;
    xprs = 0.0;

    double xdpdrho = 0.0;
    double xdedt = 0.0;
    double xdpdt = 0.0;
    double xdedrho = 0.0;
    int npoints =1;

    // printf("rho = %.15e , Ye = %.15e , T = %.15e\n",rho,ye,*temp);

    nuc_eos_m_kt1_dpdrhot_dedrhot(&npoints,&rho,temp,&ye,eps,&xeps2,press,
                    &xdpdrho,&xdpdt,&xdedrho,&xdedt,&prec,keyerr,&anyerr);

    if (*keyerr != 0) printf("EOS_Omni_dEp Keyerr: %d, Temp: %g Rho: %g Eps: %g Ye: %g\n",*keyerr,*temp,rho,*eps,ye);
    
    *dedrho = (xeps2 / rho) * xdedrho; //because of log
    *dedt   = (xeps2 / *temp) * xdedt; 
    *dpdrho = (*press / rho) * xdpdrho;
    *dpdt   = (*press / *temp) * xdpdt;

    // printf("%.15e %.15e %.15e %.15e %.15e %.15e\n",*eps,*press,*dpdrho,*dpdt,*dedrho,*dedt);

}
