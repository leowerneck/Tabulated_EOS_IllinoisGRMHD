
#ifndef NUC_EOS_HH
#define NUC_EOS_HH

//Unit changes from cgs (and temp in MeV) to HARM units
#define HAVEGR 1

#ifndef MAX
# define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

//#define MASS_BH 3.0

#define NTABLES 19
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define TEMPGF 8.617330337217213e-11
#define TEMPCONV 1.0
#define INVLENGTH 147651.770773
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38
#define INVTEMPGF  1.1604522060401007e10

#define Q_CONV 1.42134688491e-56// conversion from eV/(cm3s) to code units (for neutrinos) Assumes 1 solar mass BH (set by LENGTHGF)
#define R_CONV 7.97315453692e-24// conversion from mass/(cm3s) to code units (for neutrinos) Assumes 1 solar mass BH (set by LENGTHGF)

#ifndef GAMMA 
#define GAMMA             (gam)
#endif

#define KBERG  1.38064852e-16
#define AMU    1.6605402e-24

double temp0, temp1;
double energy_shift;

// min and max values

double eos_rhomax, eos_rhomin;
double eos_tempmin, eos_tempmax;
double eos_yemin, eos_yemax;
double eos_epsmin, eos_epsmax;
double eos_pressmin, eos_pressmax;

double c2p_tempmin;
double c2p_tempmax;

// table key
// 0 logpress 
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma
enum eos_var {i_logpress=0, i_logenergy, i_entropy, i_munu, i_cs2, i_dedt,
                i_dpdrhoe, i_dpderho, i_muhat, i_mu_e, i_mu_p, i_mu_n, i_Xa,
                i_Xh, i_Xn, i_Xp, i_Abar, i_Zbar, i_Gamma};

// table data

int nrho;
int ntemp;
int nye;
//
double * restrict alltables;
double * restrict epstable;
double * restrict presstable;
double * restrict logrho;
double * restrict logtemp;
double dlintemp,dlintempi;
double drholintempi;
double dlintempyei;
double drholintempyei;
double * restrict yes;
double dtemp, dtempi;
double drho, drhoi;
double dye, dyei;
double drhotempi;
double drhoyei;
double dtempyei;
double drhotempyei;
void nuc_eos_C_ReadTable(char* nuceos_table_name);

#endif // NUC_EOS_HH
