#ifndef NRPy_EOS_H_
#define NRPy_EOS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

// EOS struct
typedef struct _NRPyEOS_params_ {

  // Number of points
  int nrho;
  int ntemp;
  int nye;

  // Table arrays
  double *restrict alltables;
  double *restrict epstable;
  double *restrict logrho;
  double *restrict logtemp;
  double *restrict yes;

  // Minimum and maximum values of
  // rho, Ye, and T
  double eos_rhomax , eos_rhomin;
  double eos_tempmin, eos_tempmax;
  double eos_yemin  , eos_yemax;

  // Auxiliary variables
  double energy_shift;
  double temp0, temp1;
  double dlintemp, dlintempi;
  double drholintempi;
  double dlintempyei;
  double drholintempyei;
  double dtemp, dtempi;
  double drho, drhoi;
  double dye, dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;

} NRPyEOS_params;

// Table keys
#define NRPyEOS_press_key   0
#define NRPyEOS_eps_key     1
#define NRPyEOS_entropy_key 2
#define NRPyEOS_munu_key    3
#define NRPyEOS_cs2_key     4
#define NRPyEOS_depsdT_key  5
#define NRPyEOS_dPdrho_key  6
#define NRPyEOS_dPdeps_key  7
#define NRPyEOS_muhat_key   8
#define NRPyEOS_mu_e_key    9
#define NRPyEOS_mu_p_key    10
#define NRPyEOS_mu_n_key    11
#define NRPyEOS_Xa_key      12
#define NRPyEOS_Xh_key      13
#define NRPyEOS_Xn_key      14
#define NRPyEOS_Xp_key      15
#define NRPyEOS_Abar_key    16
#define NRPyEOS_Zbar_key    17
#define NRPyEOS_Gamma_key   18
#define NRPyEOS_ntablekeys  19

// Keys for table entries

// Name of the variables. This is only used to print
// information about the keys during startup
static const char table_var_names[NRPyEOS_ntablekeys][10] = {
  "logpress","logenergy","entropy","munu","cs2","dedt",
  "dpdrhoe", "dpderho", "muhat", "mu_e", "mu_p", "mu_n",
  "Xa","Xh","Xn","Xp","Abar","Zbar","Gamma"
};

// Error handling struct
typedef struct _NRPyEOS_error_report_ {
  bool error;
  int error_key;
  char message[512];
} NRPyEOS_error_report;

// Function prototypes
#include "NRPyEOS_Tabulated_headers.h"

// Helper (inline) functions
#include "NRPyEOS_Tabulated_helpers.h"

#endif // NRPy_EOS_H_
