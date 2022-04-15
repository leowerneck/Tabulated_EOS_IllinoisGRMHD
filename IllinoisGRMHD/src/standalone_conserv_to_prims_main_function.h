

#define CCTK_REAL double
#define CCTK_EQUALS(a,b) (strcmp((a),(b))==0)
#define CCTK_Equals(a,b) (strcmp((a),(b))==0)
#define cGH int

int conserv_to_prims_debug = 0;
char verbose[100];

CCTK_REAL GAMMA_SPEED_LIMIT,rho_b_atm,tau_atm, rho_b_max, Psi6threshold;

CCTK_REAL Gamma_th, K_ppoly_tab0;
CCTK_REAL rho_ppoly_tab_in[10],Gamma_ppoly_tab_in[11];

int neos;
int update_Tmunu;

#define CCTK_THORNSTRING ""
#define CCTK_WARN_ALERT  ""
#define Symmetry "none"
#define CCTK_GFINDEX3D(IGNORE,i,j,k) ((i) + cctk_lsh[0] * ((j) + cctk_lsh[1] * (k)))
#define GetRefinementLevel(cctkGH) 0

#include <stdarg.h>
#include <string.h>
#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"
#include "harm_u2p_util.c"
#include "inlined_functions.h"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"
#include "convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij.C"

int CCTK_VInfo(const char *thorn, const char *format, ...) {
  va_list ap;
  fprintf (stdout, "INFO (NOTHORN): ");
  va_start (ap, format);
  vfprintf (stdout, format, ap);
  va_end (ap);
  fprintf (stdout, "\n");
  return 0;
}

int *cctkGH;

int main(int argc, const char *argv[]) {

  if(argc != 2) {
    fprintf(stderr,"Error: Correct usage: ./driver_conserv_to_prims [filename]\n");
    exit(1);
  }

  sprintf(verbose,"essential+iteration output");
  // We use proper C++ here, for file I/O later.
  using namespace std;

  ifstream myfile;
  char filename[100];
  sprintf(filename,"%s",argv[1]);
  myfile.open (filename, ios::in | ios::binary);
  if(myfile.fail()) {
    fprintf(stderr,"Error: file %s cannot be opened.\n",filename);
    exit(1);
  }
  //myfile.open ("data.bin", ios::out | ios::binary);
  int cctk_lsh[3];
  myfile.read((char*)cctk_lsh, 3*sizeof(int));

  myfile.read((char*)&GAMMA_SPEED_LIMIT, 1*sizeof(CCTK_REAL));

  myfile.read((char*)&rho_b_max, 1*sizeof(CCTK_REAL));
  myfile.read((char*)&rho_b_atm, 1*sizeof(CCTK_REAL));
  myfile.read((char*)&tau_atm, 1*sizeof(CCTK_REAL));

  myfile.read((char*)&Psi6threshold, 1*sizeof(CCTK_REAL));

  myfile.read((char*)&update_Tmunu, 1*sizeof(int));

  myfile.read((char*)&neos,     1*sizeof(int));

  if(neos > 11 || neos < 1) {
    fprintf(stderr,"ERROR: neos = %d too large or too small\n",neos);
    exit(1);
  }

  myfile.read((char*)&Gamma_th,               1*sizeof(CCTK_REAL));
  myfile.read((char*)&K_ppoly_tab0,           1*sizeof(CCTK_REAL));
  myfile.read((char*)Gamma_ppoly_tab_in,   neos*sizeof(CCTK_REAL));
  myfile.read((char*)rho_ppoly_tab_in, (neos-1)*sizeof(CCTK_REAL));

  int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  CCTK_REAL *x = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *y = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *z = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)x,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)y,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)z,   (fullsize)*sizeof(CCTK_REAL));

  // Should probably output these:
  CCTK_REAL *failure_checker = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtt = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTty = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTyy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTyz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTzz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *alp = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gyy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gyz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gzz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *psi_bssn = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char *)failure_checker, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTtt, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTtx, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTty, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTtz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTxx, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTxy, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTxz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTyy, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTyz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)eTzz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)alp, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gxx, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gxy, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gxz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gyy, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gyz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)gzz, fullsize*sizeof(CCTK_REAL));
  myfile.read((char *)psi_bssn, fullsize*sizeof(CCTK_REAL));



  CCTK_REAL *phi_bssn = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtyy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtyz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtzz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)phi_bssn, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtxx, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtxy, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtxz, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtyy, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtyz, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtzz, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *gtupxx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupxy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupxz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupyy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupyz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupzz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)gtupxx, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtupxy, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtupxz, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtupyy, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtupyz, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)gtupzz, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *betax = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *betay = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *betaz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)betax, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)betay, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)betaz, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *lapm1 = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)lapm1, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *tau = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_x = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_y = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_z = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)tau,      (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)mhd_st_x, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)mhd_st_y, (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)mhd_st_z, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *rho_star = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)rho_star, (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *Bx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *By = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *Bz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)Bx,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)By,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)Bz,   (fullsize)*sizeof(CCTK_REAL));

  // v^i = u^i / u^0 <-- NOT VALENCIA VELOCITY
  CCTK_REAL *vx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *vy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *vz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)vx,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)vy,   (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)vz,   (fullsize)*sizeof(CCTK_REAL));

  CCTK_REAL *P = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *rho_b = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  myfile.read((char*)P,    (fullsize)*sizeof(CCTK_REAL));
  myfile.read((char*)rho_b,(fullsize)*sizeof(CCTK_REAL));

  //int checker=1063;
  int checker;
  myfile.read((char*)&checker,sizeof(int));
  if(checker != 1063) {
    fprintf(stderr,"MAGIC NUMBER FAILED. DATA FILE READIN ERROR.\n");
    exit(1);
  }

  myfile.close();

  // HERE WE USE _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
  CCTK_REAL *rho_star_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_x_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_y_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_z_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *tau_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
