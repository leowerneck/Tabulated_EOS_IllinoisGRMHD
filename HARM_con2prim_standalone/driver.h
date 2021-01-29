#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
//#include "decs.h"
#include "nuc_eos.h"
//#include "EOS_functions.c"
#include "units.h"

//#include "utoprim_2d_new.c"
/* your choice of floating-point data type */
#define FTYPE double 

#define RHO      (0)
#define UU       (1)
#define u1_con   (2)
#define u2_con   (3)
#define u3_con   (4)
#define B1_con   (5)
#define B2_con   (6)
#define B3_con   (7)
#define YE       (8)
#define TEMP     (9)
#define EPS      (10)
#define PRESS    (11)
#define WLORENTZ (12)

/* for conserved variables */
#define QCOV0   1
#define QCOV1   2
#define QCOV2   3
#define QCOV3   4

#define UTCON1  (u1_con)
#define UTCON2  (u2_con)
#define UTCON3  (u3_con)

#define U1  (u1_con)
#define U2  (u2_con)
#define U3  (u3_con)

#define BCON1   (B1_con)
#define BCON2   (B2_con)
#define BCON3   (B3_con)

// some constants (required by EOS_Omni
// in standalone version)
extern const double kBerg;
extern const double amu;

// conversion factors between cgs and M_Sun = c = G = 1
// (temperature between K and MeV)
// see EOS_Omni/doc/units.py
extern const double rho_gf;
extern const double press_gf;
extern const double eps_gf;
extern const double temp_gf;

// Inverses of the numbers above
extern const double inv_rho_gf;
extern const double inv_press_gf;
extern const double inv_eps_gf;
extern const double inv_temp_gf;


struct metric {
  double lo[4][4];
  double up[4][4];
  double lo_det;
  double alpha;
  double up_det;
  double lo_sqrt_det;
};

struct of_geom {
  double gcon[4][4] ;
  double gcov[4][4] ;
  double g ;
  double g_inv;
  double alpha;
  double beta[3];
  double ncon[4];
} ;


struct report {
  bool failed;
  bool adjust_cons;
  char err_msg[200];
  int count;
  bool retry;
  int keyerr;
};


extern void nuc_eos_m_kt0_press(const int *restrict n_in,
                                const double *restrict rho,
                                double *restrict temp,
                                const double *restrict ye,
                                const double *restrict eps,
                                double *restrict prs,
                                const double *restrict prec,
                                int *restrict keyerr,
                                int *restrict anyerr);

extern int Utoprim_2d_eos(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);

extern int Utoprim_3d_eos(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);

extern int Utoprim_2d_eos_backup(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);

extern int Utoprim_2d_eos_safe_guess(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);

extern void nuc_eos_P_from_Enthalpy(double rho, double *temp, double ye,
                                    double *enr, double *enr2, double *press,
                                    int keytemp,
                                    int *keyerr,double rfeps); 

//#define FTYPE double 
//static int Utoprim_2d_fast(FTYPE *U, struct of_geom *geom, FTYPE *prim, FTYPE *gamma, FTYPE *bsq);


//void nuc_eos_m_kt1_dpdrhoe_dpderho(const int *restrict n_in,
//                                     const double *restrict rho,
//                                     double *restrict temp,
//                                     const double *restrict ye,
//                                     const double *restrict eps,
//                                     double *restrict dpdrhoe,
//                                     double *restrict dpderho,
//                                     const double *restrict prec,
//                                     int *restrict keyerr,
//                                    int *restrict anyerr);

//void nuc_eos_m_kt0_press_cs2(const int *restrict n_in,
//                         const double *restrict rho,
//                         double *restrict temp,
//                         const double *restrict ye,
//                         const double *restrict eps,
//                        double *restrict prs,
//                         double *restrict cs2,
//                         const double *restrict prec,
//                        int *restrict keyerr,
//                        int *restrict anyerr);


//void EOS_press(double rho, double * eps, 
//        double * temp, double ye, double * prs);

//void EOS_Omni_press(double prec, double rho, double * eps, 
//    double * temp, double ye,  double * press);



//void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, 
//          double * temp_guess, double ye, double * eps, double * press, double * dPdrho, 
//          double * dPdeps, double * entr,const double prec);

//void EOS_Omni_dpdrho_dpdeps_hinv(const double rho, const double enth, double * temp, 
//       const double ye, double * eps, double * press, double * dpdrho, double * dpdeps, double * entr, 
//      int * keyerr, int stepsize);

extern void EOS_press_with_T(double rho, double * eps, 
                             double * temp, double ye, double * prs);

extern void EOS_Omni_press_with_T(double prec, double rho, double * eps, 
                                  double * temp, double ye,  double * press);

extern void EOS_press(double rho, double * eps, 
                      double * temp, double ye, double * prs);

extern void EOS_Omni_press(double prec, double rho, double * eps, 
                           double * temp, double ye,  double * press);

extern void EOS_press_cs2(double rho, double * eps, 
                          double * temp, double ye, double * prs, double * cs2);

extern void EOS_Omni_press_cs2(double prec, double rho, double * eps, 
                               double * temp, double ye,  double * press, double * cs2);


extern void EOS_P_from_hrho_dPdrho_dPdeps(const double rho, const double enth, 
                                          double * temp_guess, double ye, double * eps, double * press, double * dPdrho, 
                                          double * dPdeps, double * entr,const double prec);

extern void EOS_Omni_dpdrho_dpdeps_hinv(const double rho, const double enth, double * temp, 
                                        const double ye, double * eps, double * press, double * dpdrho, double * dpdeps, double * entr, 
                                        int * keyerr, int stepsize);


extern void EOS_press_testing_eos(int eoskey, int keytemp, double rho, double * eps, 
                                  double * temp, double ye, double * prs, int * keyerr);

extern void EOS_Omni_press_testing_eos(int eoskey, int keytemp, double prec, double rho, double * eps, 
                                       double * temp, double ye, double * press, int * keyerr);


extern void nuc_eos_m_kt0_press(const int *restrict n_in,
                                const double *restrict rho,
                                double *restrict temp,
                                const double *restrict ye,
                                const double *restrict eps,
                                double *restrict prs,
                                const double *restrict prec,
                                int *restrict keyerr,
                                int *restrict anyerr);



extern void nuc_eos_P_from_Enthalpy(double rho, double *temp, double ye,
                                    double *enr, double *enr2, double *press,
                                    int keytemp,
                                    int *keyerr,double rfeps); 




extern  void nuc_eos_m_kt1_dpdrhoe_dpderho(const int *restrict n_in,
                                           const double *restrict rho,
                                           double *restrict temp,
                                           const double *restrict ye,
                                           const double *restrict eps,
                                           double *restrict dpdrhoe,
                                           double *restrict dpderho,
                                           const double *restrict prec,
                                           int *restrict keyerr,
                                           int *restrict anyerr);

extern  void nuc_eos_m_kt0_press_cs2(const int *restrict n_in,
                                     const double *restrict rho,
                                     double *restrict temp,
                                     const double *restrict ye,
                                     const double *restrict eps,
                                     double *restrict prs,
                                     double *restrict cs2,
                                     const double *restrict prec,
                                     int *restrict keyerr,
                                     int *restrict anyerr);


extern void EOS_Omni_dpdrho_dpdt_dedrho_dedt(int eoskey, double rho, double * temp, 
                                             double ye, double * eps, double * press, double * dpdrho, double * dpdt, double * dedrho, 
                                             double * dedt, int * keyerr);

extern void EOS_EP_dEdr_dEdt_dPdr_dPdt(double * x, const double D, 
                                       double * Eprim, double * Pprim, double * dEdrho, double * dEdt, double * dPdrho, 
                                       double * dPdt, double ye, int eoskey);

extern void EP_dEdgamma_dEdW_dEdT(double * Eprim, double * Pprim, double * dEdgamma, double * dEdW, double * dEdT, 
                                  double * dpdrho, double * dpdT, double * x, const double D, double ye, int eoskey);

