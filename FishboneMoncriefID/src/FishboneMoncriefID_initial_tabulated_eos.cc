#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "FishboneMoncriefID.h"

namespace nuc_eos_private {
  extern int nrho;
  extern CCTK_REAL *restrict logrho;
}

// Find interpolation index using Bisection root-finding algorithm:
static inline int bisect_left(const CCTK_REAL rrbar, const int numlines_in_file, const CCTK_REAL *restrict rbar_arr) {

  int x1 = 0;
  int x2 = numlines_in_file-1;
  CCTK_REAL y1 = rrbar-rbar_arr[x1];
  CCTK_REAL y2 = rrbar-rbar_arr[x2];
  if(y1*y2 >= 0) {
    CCTK_VINFO("INTERPOLATION BRACKETING ERROR INSIDE FUNCTION: %s",__func__);
    CCTK_VINFO("Indices : %d %d",x1,x2);
    CCTK_VINFO("Values  : %e %e",exp(rbar_arr[x1]),exp(rbar_arr[x2]));
    CCTK_VINFO("Interval: %e %e",exp(y1),exp(y2));
    CCTK_VINFO("Point   : %e",exp(rrbar));
    CCTK_ERROR("ABORTING");
  }
  for(int i=0;i<numlines_in_file;i++) {
    int x_midpoint = (x1+x2)/2;
    CCTK_REAL y_midpoint = rrbar-rbar_arr[x_midpoint];
    if(y_midpoint*y1 < 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if( abs(x2-x1) == 1 ) {
      // Always return the left value
      return x1;
    }
  }
  CCTK_ERROR("INTERPOLATION BRACKETING ERROR: DID NOT CONVERGE.");
}

// Reference: https://en.wikipedia.org/wiki/Linear_interpolation
inline void linear_interp_1D( const int nx,
                              const CCTK_REAL *restrict x_arr,
                              const CCTK_REAL *restrict y_arr,
                              const CCTK_REAL *restrict z_arr,
                              const CCTK_REAL x_star,
                              CCTK_REAL *restrict y_star,
                              CCTK_REAL *restrict z_star) {
  // Find the index to the left
  int idx = bisect_left(x_star,nx,x_arr);
  // Set (x0,y0,z0) and (x1,y1,z1)
  const CCTK_REAL x0 = x_arr[idx  ];
  const CCTK_REAL y0 = y_arr[idx  ];
  const CCTK_REAL z0 = z_arr[idx  ];
  const CCTK_REAL x1 = x_arr[idx+1];
  const CCTK_REAL y1 = y_arr[idx+1];
  const CCTK_REAL z1 = z_arr[idx+1];
  // Perform the interpolation. Note that:
  //
  // y_star = y0 + (x_star-x0)*(y1-y0)/(x1-x0)
  //        = y0 - y0*(x_star-x0)/(x1-x0) + y1*(x_star-x0)/(x1-x0)
  //    .--------------------------------.
  // => | y_star = y0*(1 - aux) + y1*aux | ,
  //    .--------------------------------.
  // where
  //    .---------------------------.
  //    | aux = (x_star-x0)/(x1-x0) | .
  //    .---------------------------.
  const CCTK_REAL aux = (x_star-x0)/(x1-x0);
  *y_star = y0*(1.0-aux) + y1*aux;
  *z_star = z0*(1.0-aux) + z1*aux;
}

extern "C"
void FishboneMoncriefID_initial_tabulated_eos(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // First compute maximum pressure and density
  CCTK_REAL h_max;
  {
    CCTK_REAL hm1;
    CCTK_REAL xcoord = r_at_max_density;
    CCTK_REAL ycoord = 0.0;
    CCTK_REAL zcoord = 0.0;
    {
#include "FishboneMoncriefID_hm1.h"
    }
    h_max = hm1+1;
  }

  // Step 1: Compute atmosphere pressure and density
  CCTK_REAL atmosphere_P,atmosphere_eps,atmosphere_S,atmosphere_h;
  const int npoints   = 1;
  const int eos_key   = EOS_Omni_GetHandle("nuc_eos");
  CCTK_REAL precision = 1e-10;
  CCTK_REAL T_in = atmosphere_temperature;
  if( initialize_entropy ) {
    int key_err    = 0;
    int any_err    = 0;
    int temp_key   = 1; // 1 means "have temperature"
    CCTK_REAL __attribute__((unused)) cs2, dedt, dpderho, dpdrhoe, munu;
    EOS_Omni_short(eos_key,
                   temp_key,
                   precision,
                   npoints,
                   &atmosphere_rho,
                   &atmosphere_eps,
                   &T_in,
                   &atmosphere_Y_e,
                   &atmosphere_P,
                   &atmosphere_S,
                   &cs2,&dedt,&dpderho,&dpdrhoe,&munu,
                   &key_err,&any_err);
    if( any_err ) CCTK_VWARN(CCTK_WARN_ALERT,"Problem in EOS_Omni_short when setting atmosphere quantities. key_err = %d",key_err);
  }
  else {
    int key_err    = 0;
    int any_err    = 0;
    int temp_key   = 1; // 1 means "have temperature"
    EOS_Omni_press(eos_key,
                   temp_key,
                   precision,
                   npoints,
                   &atmosphere_rho,
                   &atmosphere_eps,
                   &T_in,
                   &atmosphere_Y_e,
                   &atmosphere_P,
                   &key_err,&any_err);
    if( any_err ) CCTK_VWARN(CCTK_WARN_ALERT,"Problem in EOS_Omni_press when setting atmosphere quantities. key_err = %d",key_err);
  }

  atmosphere_h = 1 + atmosphere_eps + atmosphere_P/atmosphere_rho;

  if( verbosity_level > 0 ) {
    CCTK_VINFO("Atmosphere parameters:");
    CCTK_VINFO("  - Density          : %e",atmosphere_rho);
    CCTK_VINFO("  - Electron fraction: %e",atmosphere_Y_e);
    CCTK_VINFO("  - Temperature      : %e (=%e)",atmosphere_temperature,T_in);
    CCTK_VINFO("  - Pressure         : %e",atmosphere_P);
    CCTK_VINFO("  - Energy           : %e",atmosphere_eps);
    CCTK_VINFO("  - Entropy          : %e",atmosphere_S);
    CCTK_VINFO("  - Enthalpy         : %e",atmosphere_h);
  }

  // Step 2: Slice the table at constant Y_e and S
  // Step 2.a: Allocate memory for log(rho), log(T), and log(h)
  CCTK_REAL *lr_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*nuc_eos_private::nrho);
  CCTK_REAL *lt_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*nuc_eos_private::nrho);
  CCTK_REAL *lh_arr = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*nuc_eos_private::nrho);

  int max_idx = 0;
  CCTK_REAL hm1_min = 0;
  // Step 2.b: Loop over table densities computing T from (rho,Ye,S); compute h
  for(int ir=0;ir<nuc_eos_private::nrho;ir++) {

    // Step 3.a: Set the local density
    const CCTK_REAL lr   = nuc_eos_private::logrho[ir];
    const CCTK_REAL rhoL = exp(lr);

    // Step 3.c: Given (rho,Y_e,S), determine T, P, and eps
    int key_err            = 0;
    int any_err            = 0;
    int temp_key           = 2; // 2 means "compute T from S"
    CCTK_REAL Y_eL         = disk_Y_e;
    CCTK_REAL entropyL     = disk_entropy;
    CCTK_REAL temperatureL = 1.0; // Temperature guess
    CCTK_REAL pressL, epsL;
    CCTK_REAL __attribute__((unused)) cs2, dedt, dpderho, dpdrhoe, munu;
    EOS_Omni_short(eos_key,
                   temp_key,
                   precision,
                   npoints,
                   &rhoL,&epsL,&temperatureL,&Y_eL,&pressL,&entropyL,
                   &cs2,&dedt,&dpderho,&dpdrhoe,&munu,
                   &key_err,&any_err);
    if( any_err ) CCTK_VWARN(CCTK_WARN_ALERT,"Problem in EOS_Omni_short when slicing the table. key_err = %d",key_err);

    // Step 3.d: Compute the local enthalpy
    const CCTK_REAL h = 1 + epsL + pressL/rhoL;

    // Step 3.e: Update the rho, T, h arrays
    lr_arr[ir] = lr;
    lt_arr[ir] = log(temperatureL);
    lh_arr[ir] = log(h);
    if( verbosity_level > 1 ) CCTK_VINFO("Slice: %d (%d) %.15e %.15e %.15e",max_idx,nuc_eos_private::nrho,rhoL,temperatureL,h);
    max_idx++;
    if(  ir == 0  ) hm1_min = h-1;
    if( h > h_max ) break;
  }

  if( verbosity_level > 0 ) {
    CCTK_VINFO("Minimum enthalpy in disk: %.15e",hm1_min+1);
    CCTK_VINFO("Maximum enthalpy in disk: %.15e",h_max);
  }

  //#pragma omp parallel for
  for(CCTK_INT k=0;k<cctk_lsh[2];k++) {
    for(CCTK_INT j=0;j<cctk_lsh[1];j++) {
      for(CCTK_INT i=0;i<cctk_lsh[0];i++) {

        // Step 4: Set local gridfunction index
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        // Step 5: Set local coordinates (x,y,z,r)
        CCTK_REAL xcoord = x[idx];
        CCTK_REAL ycoord = y[idx];
        CCTK_REAL zcoord = z[idx];
        CCTK_REAL rr     = r[idx];

        // Step 6: Compute metric quantities
        FishboneMoncriefID_KerrSchild(cctkGH,cctk_lsh,
                                      i,j,k,
                                      x,y,z,
                                      alp,betax,betay,betaz,
                                      gxx,gxy,gxz,gyy,gyz,gzz,
                                      kxx,kxy,kxz,kyy,kyz,kzz);

        // Step 7: Compute Hydro quantities
        CCTK_REAL hm1;
        bool set_to_atmosphere=false;
        if(rr > r_in) {
          // Step 7.a: Compute enthalpy
          {
#include "FishboneMoncriefID_hm1.h"
          }
          if(hm1 > hm1_min) {
            // Step 7.b: Given h, interpolate rho and T
            CCTK_REAL lh = log(hm1+1);
            CCTK_REAL lr,lt;
            linear_interp_1D(max_idx,lh_arr,lr_arr,lt_arr,lh,&lr,&lt);

            // Step 7.c: Compute P and eps
            int key_err   = 0;
            int any_err   = 0;
            int temp_key  = 1; // 1 means "have temperature"
            const CCTK_REAL rhoL   = exp(lr);
            const CCTK_REAL Y_eL   = disk_Y_e;
            CCTK_REAL temperatureL = exp(lt);
            CCTK_REAL pressL,epsL;
            EOS_Omni_press(eos_key,
                           temp_key,
                           precision,
                           npoints,
                           &rhoL,&epsL,&temperatureL,&Y_eL,&pressL,
                           &key_err,&any_err);

            // Step 7.d: Write to hydro gridfunctions
            rho        [idx] = rhoL;
            Y_e        [idx] = Y_eL;
            temperature[idx] = temperatureL;
            press      [idx] = pressL;
            eps        [idx] = epsL;
            if( initialize_entropy ) {
              entropy  [idx] = disk_entropy;
            }

            // Step 7.e: Compute velocities
            FishboneMoncriefID_velocities(cctkGH,cctk_lsh,
                                          i,j,k,
                                          x,y,z,
                                          velx,vely,velz);
          } else {
            set_to_atmosphere=true;
          }
        } else {
          set_to_atmosphere=true;
        }
        // 7.f: Outside the disk? Set to atmosphere all hydrodynamic variables!
        if(set_to_atmosphere) {
          rho        [idx] = atmosphere_rho;
          Y_e        [idx] = atmosphere_Y_e;
          temperature[idx] = atmosphere_temperature;
          press      [idx] = atmosphere_P;
          eps        [idx] = atmosphere_eps;
          w_lorentz  [idx] = 1.0;
          velx       [idx] = 0.0;
          vely       [idx] = 0.0;
          velz       [idx] = 0.0;
          if( initialize_entropy ) {
            entropy  [idx] = atmosphere_S;
          }
        }
      }
    }
  }

  // Step 8: Free memory
  free(lr_arr);
  free(lt_arr);
  free(lh_arr);
}
