// Rewrite of Tmunu.F90 in C++

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//------------------------------------------
//---------- IllinoisGRMHD stuff -----------
//------------------------------------------
#include "../IllinoisGRMHD_headers.h"

//------------------------------------------
//------------- EOS_Omni stuff -------------
//------------------------------------------
#define INVRHOGF (6.17714470405638e17)
#define PRESSGF  (1.80123683248503e-39)
#define RHOGF    (1.61887093132742e-18)
namespace nuc_eos {
  extern CCTK_REAL eos_rhomax;
}

//------------------------------------------
//----------- ZelmaniLeak stuff ------------
//------------------------------------------
#define PRESS_NU_CONSTANT (3.52127727)
#define PI2               (9.86960440108936)
#define PI4               (97.4090910340024)

namespace ZelmaniLeak {

  void ZelmaniLeak_Tmunu( const igm_eos_parameters eos, CCTK_ARGUMENTS ) {

    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // Find the maximum density in the grid
    CCTK_INT  npoints      = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    CCTK_REAL rho_max_grid = rho_b[0];
    for( int i=1;i<npoints;i++ ) {
      CCTK_REAL rho_aux = rho_b[i];
      if( rho_aux > rho_max_grid ) rho_max_grid = rho_aux;
    }

    // Maximum allowed value of rho
    const CCTK_REAL rho_max_cgs = rho_max_grid*INVRHOGF;

    // If below
    if( rho_max_cgs < pnu_rho_start/10.0 ) return;

    const CCTK_REAL pnuconst = PRESS_NU_CONSTANT * PRESSGF;
    const CCTK_REAL F3const  = 7.0*PI4/60.0;

    CCTK_REAL enufac;
    if( include_enu_in_tmunu != 0 )
      enufac = 1.0;
    else
      enufac = 0.0;

    const CCTK_INT imin = 0, imax = cctk_lsh[0];
    const CCTK_INT jmin = 0, jmax = cctk_lsh[1];
    const CCTK_INT kmin = 0, kmax = cctk_lsh[2];

#pragma omp parallel for
    for(int k=kmin;k<kmax;k++) {
      for(int j=jmin;j<jmax;j++) {
        for(int i=imin;i<imax;i++) {
          const CCTK_INT index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          //---------------------------------------------------------
          //------------------ IllinoisGRMHD stuff ------------------
          //---------------------------------------------------------
          CCTK_REAL METRIC[NUMVARS_FOR_METRIC];
          METRIC[PHI   ] = phi_bssn[index];
          METRIC[GXX   ] = gtxx    [index];
          METRIC[GXY   ] = gtxy    [index];
          METRIC[GXZ   ] = gtxz    [index];
          METRIC[GYY   ] = gtyy    [index];
          METRIC[GYZ   ] = gtyz    [index];
          METRIC[GZZ   ] = gtzz    [index];
          METRIC[LAPM1 ] = lapm1   [index];
          METRIC[SHIFTX] = betax   [index];
          METRIC[SHIFTY] = betay   [index];
          METRIC[SHIFTZ] = betaz   [index];
          METRIC[GUPXX ] = gtupxx  [index];
          METRIC[GUPYY ] = gtupyy  [index];
          METRIC[GUPZZ ] = gtupzz  [index];
          METRIC[GUPXY ] = gtupxy  [index];
          METRIC[GUPXZ ] = gtupxz  [index];
          METRIC[GUPYZ ] = gtupyz  [index];

          CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
          SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);

          CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
          METRIC_PHYS[GXX  ] = METRIC[GXX  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GXY  ] = METRIC[GXY  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GXZ  ] = METRIC[GXZ  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GYY  ] = METRIC[GYY  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GYZ  ] = METRIC[GYZ  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GZZ  ] = METRIC[GZZ  ]*METRIC_LAP_PSI4[PSI4 ];
          METRIC_PHYS[GUPXX] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPXY] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPXZ] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPYY] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPYZ] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4];
          METRIC_PHYS[GUPZZ] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4];
          //---------------------------------------------------------

          // Compute munu from rho, Ye, T
          CCTK_REAL xrho  = rho_b[index];
          CCTK_REAL xye   = igm_Ye[index];
          CCTK_REAL xtemp = igm_temperature[index];
          CCTK_REAL xmunu = 0.0;
          get_munu_from_rho_Ye_and_T( eos,xrho,xye,xtemp, &xmunu );

          // Continue with ZelmaniLeak stuff
          CCTK_REAL eta          = xmunu/xtemp;
          CCTK_REAL half_eta_sqr = 0.5*eta*eta;
          CCTK_REAL F3           = F3const + half_eta_sqr*(PI2 + half_eta_sqr);
          CCTK_REAL pnu          = F3 * pnuconst * pow(xtemp,4.0) * exp( -pnu_rho_start*RHOGF/xrho );

          // Now first compute the Valencia velocity
          CCTK_REAL velx  = (vx[index] + METRIC[SHIFTX])*METRIC_LAP_PSI4[LAPSEINV];
          CCTK_REAL vely  = (vy[index] + METRIC[SHIFTY])*METRIC_LAP_PSI4[LAPSEINV];
          CCTK_REAL velz  = (vz[index] + METRIC[SHIFTZ])*METRIC_LAP_PSI4[LAPSEINV];

          // Then lower the index
          CCTK_REAL vel_x = METRIC_PHYS[GXX]*velx + METRIC_PHYS[GXY]*vely + METRIC_PHYS[GXZ]*velz;
          CCTK_REAL vel_y = METRIC_PHYS[GXY]*velx + METRIC_PHYS[GYY]*vely + METRIC_PHYS[GYZ]*velz;
          CCTK_REAL vel_z = METRIC_PHYS[GXZ]*velx + METRIC_PHYS[GYZ]*vely + METRIC_PHYS[GZZ]*velz;

          // We must also compute beta_{i}
          CCTK_REAL beta_x = METRIC_PHYS[GXX]*METRIC[SHIFTX] + METRIC_PHYS[GXY]*METRIC[SHIFTY] + METRIC_PHYS[GXZ]*METRIC[SHIFTZ];
          CCTK_REAL beta_y = METRIC_PHYS[GXY]*METRIC[SHIFTX] + METRIC_PHYS[GYY]*METRIC[SHIFTY] + METRIC_PHYS[GYZ]*METRIC[SHIFTZ];
          CCTK_REAL beta_z = METRIC_PHYS[GXZ]*METRIC[SHIFTX] + METRIC_PHYS[GYZ]*METRIC[SHIFTY] + METRIC_PHYS[GZZ]*METRIC[SHIFTZ];

          // Now compute the Lorentz factor
          CCTK_REAL vsq = velx*vel_x + vely*vel_y + velz*vel_z;
          CCTK_REAL W   = 1.0/sqrt(1.0 - vsq);

          // Next compute beta^{2} := beta^{i}beta_{i}
          CCTK_REAL beta2 = METRIC[SHIFTX]*beta_x + METRIC[SHIFTY]*beta_y + METRIC[SHIFTZ]*beta_z;

          CCTK_REAL tmunutemp = W*W*(3.0*pnu*enufac + pnu);

          CCTK_REAL u_t = -METRIC_LAP_PSI4[LAPSE] + velx*beta_x + vely*beta_y + velz*beta_z;
          CCTK_REAL u_x = vel_x;
          CCTK_REAL u_y = vel_y;
          CCTK_REAL u_z = vel_z;

          // Finally, add to T^{\mu\nu}
          eTtt[index] += tmunutemp * u_t * u_t + pnu*( beta2 - SQR(METRIC_LAP_PSI4[LAPSE]) );

          eTtx[index] += tmunutemp * u_t * u_x + pnu*beta_x;
          eTty[index] += tmunutemp * u_t * u_y + pnu*beta_y;
          eTtz[index] += tmunutemp * u_t * u_z + pnu*beta_z;

          eTxx[index] += tmunutemp * u_x * u_x + pnu*METRIC_PHYS[GXX];
          eTyy[index] += tmunutemp * u_y * u_y + pnu*METRIC_PHYS[GYY];
          eTzz[index] += tmunutemp * u_z * u_z + pnu*METRIC_PHYS[GZZ];

          eTxy[index] += tmunutemp * u_x * u_y + pnu*METRIC_PHYS[GXY];
          eTxz[index] += tmunutemp * u_x * u_z + pnu*METRIC_PHYS[GXZ];
          eTyz[index] += tmunutemp * u_y * u_z + pnu*METRIC_PHYS[GYZ];        

        } // loop i
      } // loop j
    } // loop k
  }

} // namespace ZelmaniLeak
