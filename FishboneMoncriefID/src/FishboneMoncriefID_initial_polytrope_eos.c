#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "FishboneMoncriefID.h"

void FishboneMoncriefID_initial_polytrope_eos(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // First compute maximum pressure and density
  CCTK_REAL P_max, rho_max;
  {
    CCTK_REAL hm1;
    CCTK_REAL xcoord = r_at_max_density;
    CCTK_REAL ycoord = 0.0;
    CCTK_REAL zcoord = 0.0;
    {
#include "FishboneMoncriefID_hm1.h"
    }
    rho_max = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) );
    P_max   = kappa * pow(rho_max, gamma);
  }

  // We enforce units such that rho_max = 1.0; if these units are not obeyed, then
  //    we error out. If we did not error out, then the value of kappa used in all
  //    EOS routines would need to be changed, and generally these appear as
  //    read-only parameters.
  if(fabs(P_max/rho_max - kappa) > 1e-8) {
    CCTK_VINFO("Error: To ensure that P = kappa*rho^Gamma, where rho_max = 1.0,\n");
    CCTK_VINFO("       you must set (in your parfile) the polytropic constant kappa = P_max/rho_max = %.15e\n\n",P_max/rho_max);
    CCTK_VINFO(" Needed values for kappa, for common values of Gamma:\n");
    CCTK_VINFO(" For Gamma =4/3, use kappa=K_initial=K_poly = 4.249572342020724e-03 to ensure rho_max = 1.0\n");
    CCTK_VINFO(" For Gamma =5/3, use kappa=K_initial=K_poly = 6.799315747233158e-03 to ensure rho_max = 1.0\n");
    CCTK_VINFO(" For Gamma = 2,  use kappa=K_initial=K_poly = 8.499144684041449e-03 to ensure rho_max = 1.0\n");
    CCTK_ERROR("ABORTING");
  }

#pragma omp parallel for
  for(CCTK_INT k=0;k<cctk_lsh[2];k++) for(CCTK_INT j=0;j<cctk_lsh[1];j++) for(CCTK_INT i=0;i<cctk_lsh[0];i++) {
        CCTK_INT idx = CCTK_GFINDEX3D(cctkGH,i,j,k);

        CCTK_REAL xcoord = x[idx];
        CCTK_REAL ycoord = y[idx];
        CCTK_REAL zcoord = z[idx];
        CCTK_REAL rr = r[idx];

        FishboneMoncriefID_KerrSchild(cctkGH,cctk_lsh,
                                      i,j,k,
                                      x,y,z,
                                      alp,betax,betay,betaz,
                                      gxx,gxy,gxz,gyy,gyz,gzz,
                                      kxx,kxy,kxz,kyy,kyz,kzz);

        CCTK_REAL hm1;
        bool set_to_atmosphere=false;
        if(rr > r_in) {
          {
#include "FishboneMoncriefID_hm1.h"
          }
          if(hm1 > 0) {
            rho[idx] = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) ) / rho_max;
            press[idx] = kappa*pow(rho[idx], gamma);
            // P = (\Gamma - 1) rho epsilon
            eps[idx] = press[idx] / (rho[idx] * (gamma - 1.0));
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
        // Outside the disk? Set to atmosphere all hydrodynamic variables!
        if(set_to_atmosphere) {
          // Choose an atmosphere such that
          //   rho =       1e-5 * r^(-3/2), and
          //   P   = k rho^gamma
          // Add 1e-100 or 1e-300 to rr or rho to avoid divisions by zero.
          rho[idx] = 1e-5 * pow(rr + 1e-100,-3.0/2.0);
          press[idx] = kappa*pow(rho[idx], gamma);
          eps[idx] = press[idx] / ((rho[idx] + 1e-300) * (gamma - 1.0));
          w_lorentz[idx] = 1.0;
          velx[idx] = 0.0;
          vely[idx] = 0.0;
          velz[idx] = 0.0;
        }
      }

  CCTK_INT final_idx = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,cctk_lsh[1]-1,cctk_lsh[2]-1);
  if( verbosity_level > 0 ) {
    CCTK_VINFO("=====   OUTPUTS   =====");
    CCTK_VINFO("betai: %e %e %e \ngij: %e %e %e %e %e %e \nKij: %e %e %e %e %e %e\nalp: %e\n",betax[final_idx],betay[final_idx],betaz[final_idx],gxx[final_idx],gxy[final_idx],gxz[final_idx],gyy[final_idx],gyz[final_idx],gzz[final_idx],kxx[final_idx],kxy[final_idx],kxz[final_idx],kyy[final_idx],kyz[final_idx],kzz[final_idx],alp[final_idx]);
    CCTK_VINFO("rho: %.15e\nPressure: %.15e\nvx: %.15e\nvy: %.15e\nvz: %.15e",rho[final_idx],press[final_idx],velx[final_idx],vely[final_idx],velz[final_idx]);
  }
}
