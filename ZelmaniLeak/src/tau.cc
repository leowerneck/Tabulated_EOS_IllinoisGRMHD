#include <stdio.h>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <Symmetry.h>
#include <loopcontrol.h>
#include <time.h>
#include "tau.hh"

using namespace ZLtau;

extern "C" {
  void ZLtau_setup(CCTK_ARGUMENTS);
  void ZLtau_setup_local(CCTK_ARGUMENTS);
  void ZLtau_get_rays(CCTK_ARGUMENTS);
  void ZLtau_register(CCTK_ARGUMENTS);
  void ZLtau_set_origin_init(CCTK_ARGUMENTS);
}

double alpha_grid(double rmin, double rmax, int nzones,
		  double drmin, double prec);

void ZLtau_setup_local(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(CCTK_QueryGroupStorage(cctkGH,"ZelmaniLeak::tau3D")) {
    for(int i=0;i<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]*3;i++) {
      tau3D[i] = 0.0e0;
    }
  }

  if(CCTK_QueryGroupStorage(cctkGH,"ZelmaniLeak::zelmani_leak_account_local")) {
    for(int i=0;i<cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]*3;i++) {
      lum_local[i]      = 0.0e0;
      lum_int_local[i]  = 0.0e0;
      net_heat_local[i] = 0.0e0;
      heat_local[i]     = 0.0e0;
      eave_local[i] = 0.0e0;
    }
  }

  CCTK_Info(CCTK_THORNSTRING,"Resetting leakage variables to zero");

}

void ZLtau_register(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  
  
  // register 'evolved' variables so that we don't catch nans
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_int_local[0]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_int_local[1]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_int_local[2]"));

  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_local[0]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_local[1]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::lum_local[2]"));

  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::net_heat_local[0]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::net_heat_local[1]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::net_heat_local[2]"));

  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::eave_local[0]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::eave_local[1]"));
  MoLRegisterConstrained(CCTK_VarIndex("ZelmaniLeak::eave_local[2]"));

}

void ZLtau_setup(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;  

  enum eisym {full=1, bitant=2, octant=3};
  enum eisym isym;

  if(CCTK_EQUALS(symm,"octant")) {
    if(ntheta > 1)
#ifdef SYMMETRIC_OPERATORS
      *dtheta =  M_PI / (ntheta-1.5) / 2.0;
#else
      *dtheta =  M_PI / (ntheta-1) / 2.0;
#endif
    else
      *dtheta = 0.0;
    if(nphi > 1) 
      *dphi =   M_PI / (nphi  -1) / 2.0;
    else
      *dphi = 0.0;
    isym = octant;
  } else if CCTK_EQUALS(symm,"bitant") {
    if(ntheta > 1)
#ifdef SYMMETRIC_OPERATORS
      *dtheta =   M_PI / (ntheta-1.5) / 2.0; // 1.5 such that 2.*1.5 matches 2*(N-1)-1
#else
      *dtheta =   M_PI / (ntheta-1) / 2.0;
#endif
    else
      *dtheta = 0.0;
    if(nphi > 1)
#ifdef SYMMETRIC_OPERATORS
      *dphi   = 2*M_PI / (nphi  -2);
#else
      *dphi   = 2*M_PI / (nphi  -1);
#endif
    else
      *dphi = 0.0;
    isym = bitant;
  } else if CCTK_EQUALS(symm,"full") {
    if(ntheta > 1) 
      *dtheta =   M_PI / (ntheta-1);
    else
      *dtheta = 0.0;
    if(nphi > 1)
#ifdef SYMMETRIC_OPERATORS
      *dphi = 2*M_PI / (nphi  -2);
#else
      *dphi = 2*M_PI / (nphi  -1);
#endif
    else
      *dphi = 0.0;
    isym = full;
  } else {
    CCTK_WARN(0,"This symmetry is unknown");
  }
  
  // inner drad
  *drad   =  rad_max / (nrad - 1);
  // outer drad is tougher, since it's nonequidistant
  double prec = 1.0e-8;
  double alpha = alpha_grid(rad_max,rad_max_outer,nrad_outer,*drad,prec);
  if(alpha < 0.0) {
    CCTK_WARN(0,"Grid setup failed!");
  }

  // set up radial coordinates
  for(int i=0;i<nrad;i++) {
    rad[i] = i*(*drad);
    //fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }
  // outer, non-equidistant part
  double dr2 = *drad;
  for(int i=nrad;i<nrad+nrad_outer;i++) {
    rad[i] = rad[i-1] + dr2;
    dr2 = dr2 * alpha;
    //  fprintf(stdout,"rad: %5d %15.6E\n",i,rad[i]);
  }

#ifdef SYMMETRIC_OPERATORS
  if(isym == full || isym == bitant) {
    for(int i=0;i<nphi/2;i++) {
      phi[i]        = i*(*dphi);
      phi[i+nphi/2] = -phi[i];
      // fprintf(stdout,"phi: %5d %15.6E\n",i,phi[i]);
    }
  } else {
    for(int i=0;i<nphi;i++) {
      phi[i]        = i*(*dphi);
    }
  }
  if(isym == full) {
    for(int i=0;i<ntheta/2;i++) {
      theta[ntheta/2-1-i] = -(2*i+1)*(*dtheta)/2;
      theta[ntheta/2  +i] = +(2*i+1)*(*dtheta)/2;
      //fprintf(stdout,"theta: %5d %15.6E\n",i,theta[i]);
    }
  } else {
    for(int i=0;i<ntheta;i++) {
      theta[i] = (2*i-1)*(*dtheta)/2;
      //fprintf(stdout,"theta: %5d %15.6E\n",i,theta[i]);
    }
  }
#else
  for(int i=0;i<nphi;i++) {
    phi[i] = i*(*dphi);
    // fprintf(stdout,"phi: %5d %15.6E\n",i,phi[i]);
  }
  for(int i=0;i<ntheta;i++) {
    theta[i] = i*(*dtheta);
    //fprintf(stdout,"theta: %5d %15.6E\n",i,theta[i]);
  }
#endif


  // Calculate coordinates
#ifdef SYMMETRIC_OPERATORS
  if(isym == full) {
    #pragma omp parallel for
    for(int ii=0;ii<nrad+nrad_outer;ii++)
      for(int jj=0;jj<ntheta;jj++)
        for(int kk=0;kk<nphi/2;kk++)
    {
      int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
      CCTK_REAL const xrad   = rad[ii];
      CCTK_REAL const xtheta = theta[jj];
      CCTK_REAL const xphi   = phi[kk];
      CCTK_REAL xx, yy, zz;
      CCTK_REAL const xx0 = 0., yy0 = 0., zz0 = 0.;

      spher2cart (xrad,xtheta,xphi, &xx,&yy,&zz, &xx0, &yy0, &zz0);
      //    fprintf(stdout,"%5d %5d %5d %15.6E %15.6E %15.6E\n",
      //	    ii,jj,kk,xx,yy,zz);
      zi_x[iind3d] = *x0 + xx;
      zi_y[iind3d] = *y0 + yy;
      zi_z[iind3d] = *z0 + zz;

      int const iind3ds = ii+(nrad+nrad_outer)*((ntheta-1-jj)+ntheta*((kk+nphi/2) % nphi));
      zi_x[iind3ds] = *x0 - xx;
      zi_y[iind3ds] = *y0 - yy;
      zi_z[iind3ds] = *z0 - zz;
    }
  } else if(isym == bitant) {
    assert(0 && "bitant mode not yet supported"); // not sure how to handle x<0 half of domain
  } else {
    #pragma omp parallel for
    for(int ii=0;ii<nrad+nrad_outer;ii++)
      for(int jj=0;jj<ntheta;jj++)
        for(int kk=0;kk<nphi;kk++)
    {
      int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
      CCTK_REAL const xrad   = rad[ii];
      CCTK_REAL const xtheta = theta[jj];
      CCTK_REAL const xphi   = phi[kk];
      CCTK_REAL xx, yy, zz;
      CCTK_REAL const xx0 = 0., yy0 = 0., zz0 = 0.;

      spher2cart (xrad,xtheta,xphi, &xx,&yy,&zz, &xx0, &yy0, &zz0);
      //    fprintf(stdout,"%5d %5d %5d %15.6E %15.6E %15.6E\n",
      //	    ii,jj,kk,xx,yy,zz);
      zi_x[iind3d] = *x0 + xx;
      zi_y[iind3d] = *y0 + yy;
      zi_z[iind3d] = *z0 + zz;
    }
  }

  for(int ii=0;ii<nrad+nrad_outer;ii++)
    for(int jj=0;jj<ntheta;jj++)
      for(int kk=0;kk<nphi;kk++)
  {
    int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
    CCTK_REAL const xrad   = rad[ii];
    CCTK_REAL const xtheta = theta[jj];
    CCTK_REAL const xphi   = fabs(phi[kk]);
    CCTK_REAL xx, yy, zz;
    CCTK_REAL const xx0 = 0., yy0 = 0., zz0 = 0.;
    
    spher2cart (xrad,xtheta,xphi, &xx,&yy,&zz, &xx0, &yy0, &zz0);
    if(kk < nphi/2 || isym == octant) {
      xx = *x0 + xx;
      yy = *y0 + yy;
    } else {
      xx = *x0 - xx;
      yy = *y0 - yy;
    }
    if(fabs(zi_x[iind3d] - xx) > 1e-8*fabs(zi_x[iind3d] + xx) || fabs(zi_y[iind3d] - yy) > 1e-8*fabs(zi_y[iind3d] + yy) || fabs(zi_z[iind3d] - zz) > 1e-8*fabs(zi_z[iind3d] + zz)) {
      fprintf(stderr, "Unexpected coordiante (%.18e,%.18e,%.18e) instead of (%.18e,%.18e,%.18e) for sphere coords (%.18e,%.18e,%.18e) index (%d,%d,%d) of (%d,%d)\n",
      xx,yy,zz,zi_x[iind3d],zi_y[iind3d],zi_z[iind3d],xrad,xtheta,xphi,ii,jj,kk,ntheta,nphi);
    }

  } //LC_ENDLOOP3 (ZLtau_get_rays);
#else
  //  LC_LOOP3 (ZLtau_get_rays,
  //          ii,jj,kk, 0,0,0, nrad,ntheta,nphi, nrad,ntheta,nphi)
  #pragma omp parallel for
  for(int ii=0;ii<nrad+nrad_outer;ii++)
    for(int jj=0;jj<ntheta;jj++)
      for(int kk=0;kk<nphi;kk++)
  {
    int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
    CCTK_REAL const xrad   = rad[ii];
    CCTK_REAL const xtheta = jj * (*dtheta);
    CCTK_REAL const xphi   = kk * (*dphi);
    CCTK_REAL xx, yy, zz;
    
    spher2cart (xrad,xtheta,xphi, &xx,&yy,&zz, x0, y0, z0);
    //    fprintf(stdout,"%5d %5d %5d %15.6E %15.6E %15.6E\n",
    //	    ii,jj,kk,xx,yy,zz);
    zi_x[iind3d] = xx;
    zi_y[iind3d] = yy;
    zi_z[iind3d] = zz;
  } //LC_ENDLOOP3 (ZLtau_get_rays);
      //} LC_ENDLOOP3 (ZLtau_get_rays);
#endif

  // some informative output
  CCTK_VInfo (CCTK_THORNSTRING,"Radii:");
  for(int ii=0;ii<nrad+nrad_outer;ii++) {
    if(ii<1)
      fprintf(stdout,"%4d %15.6E M\n",ii,rad[ii]);
    else
      fprintf(stdout,"%4d %15.6E M dr= %15.6E\n",ii,rad[ii],rad[ii]-rad[ii-1]);
  }
  CCTK_VInfo (CCTK_THORNSTRING,"Theta angles:");
  for(int jj=0;jj<ntheta;jj++) {
    fprintf(stdout,"%4d %15.6E Pi\n",jj,theta[jj]/M_PI);
  }
  CCTK_VInfo (CCTK_THORNSTRING,"Phi angles:");
  for(int kk=0;kk<nphi;kk++) {
    fprintf(stdout,"%4d %15.6E Pi\n",kk,phi[kk]/M_PI);
  }

  *have_interp_data = 0;

}

void ZLtau_get_rays(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // don't do anything if we are not interested in doing anything
  if(!force_tau_interp) 
    if(!do_tau || !(*in_prebounce || *bounce)) return;

  if( (cctk_iteration-1) % update_tau_every != 0) return;
  CCTK_Info(CCTK_THORNSTRING, "Interpolating for tau");
  // printf("it: %d\n",cctk_iteration);

  clock_t start, end;
  double elapsed;
  start = clock();

  // Interpolate
  {
    int const interp_handle = CCTK_InterpHandle (interpolator);
    assert (interp_handle >= 0);
    int const options_handle =
      Util_TableCreateFromString (interpolator_options);
    assert (options_handle >= 0);
    int const coords_handle = CCTK_CoordSystemHandle (coordinate_system);
    assert (coords_handle >= 0);
    
    // interpolate ONLY points on proc 0, then
    // do an MPI broadcast of the results
    int npoints = 0;
    if(CCTK_MyProc(cctkGH) == 0) {
      npoints = (nrad+nrad_outer)*ntheta*nphi;
    }

    void const *const interp_coords[] = { zi_x, zi_y, zi_z };
    
    CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex ("HydroBase::rho"),
      CCTK_VarIndex ("HydroBase::temperature"),
      CCTK_VarIndex ("HydroBase::Y_e"),
      CCTK_VarIndex ("ZelmaniLeak::ZLtau_grr")
    };
    int ninputs =
      sizeof input_array_indices / sizeof *input_array_indices;

    CCTK_INT const output_array_types[] = {
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
    };
    assert (sizeof output_array_types / sizeof *output_array_types == ninputs);

    void *const output_arrays[] = { 
      zi_rho, zi_temp, zi_ye, zi_grr
    };

    assert (sizeof output_arrays / sizeof *output_arrays == ninputs);

    int const ierr =
      CCTK_InterpGridArrays (cctkGH, 3,
                             interp_handle, options_handle, coords_handle,
                             npoints, CCTK_VARIABLE_REAL, interp_coords,
                             ninputs, input_array_indices,
                             ninputs, output_array_types, output_arrays);
    assert (not ierr);

    Util_TableDestroy (options_handle);
  }

  // now proc 0 broadcasts the interpolation results
  {
    int npoints = (nrad+nrad_outer)*ntheta*nphi;
    int err = MPI_Bcast((void*)zi_rho,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)zi_temp,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)zi_ye,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    err += MPI_Bcast((void*)zi_grr,npoints,MPI_DOUBLE,0,MPI_COMM_WORLD);
    assert(err==0);
  }

 // calculate the line element
#pragma omp parallel for
  for (int kk=0; kk<nphi; ++kk) {
    for (int jj=0; jj<ntheta; ++jj) {
      for (int ii=0; ii<nrad+nrad_outer-1;ii++) {
        int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
	zi_ds[iind3d] = sqrt(zi_grr[iind3d]) * (rad[ii+1]-rad[ii]);
      }
      // outer boundary; assume that metric does not change
      int const iind3d = (nrad+nrad_outer-1)+(nrad+nrad_outer)*(jj+ntheta*kk);
      int const ii = nrad+nrad_outer-1;
      zi_ds[iind3d] = sqrt(zi_grr[iind3d]) * (rad[ii]-rad[ii-1]);
    }
  }

#if 0
  if(CCTK_MyProc(cctkGH) == 0) {
    for(int ii=0;ii<nrad;ii++) {
      int kk=5;int jj=20;
      int const iind3d = ii+(nrad+nrad_outer)*(jj+ntheta*kk);
      fprintf(stderr,"%5d %15.6E %15.6E %15.6E %15.6E %15.6E\n",ii,rad[ii],zi_ds[iind3d],
      	      sqrt(zi_grr[iind3d])*(*drad),zi_grr[iind3d],pow(zi_ds[iind3d]/(*drad),2.0));
    }
  }

  abort();
#endif

  end = clock();
  elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
  CCTK_VInfo (CCTK_THORNSTRING,"time used in proc 0: %5.2f",elapsed);
  CCTK_Info(CCTK_THORNSTRING, "Done interpolating for tau");
  *have_interp_data = 1;

  return;
}

void ZLtau_set_origin_init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *x0 = 0.0;
  *y0 = 0.0;
  *z0 = 0.0;

  *heating_efficiency = 0.0;

  for(int i=0;i<(nrad+nrad_outer)*nphi*ntheta;i++) {
    zi_rho[i]  = 0.0;
    zi_eps[i]  = 0.0;
    zi_temp[i]  = 0.0;
    zi_ye[i]  = 0.0;
    zi_grr[i]  = 0.0;
    zi_ds[i]  = 0.0;
  }

  for(int i=0;i<(nrad+nrad_outer)*nphi*ntheta*3;i++) {
      zi_xiross[i] = 0.0;
      zi_tauruff[i] = 0.0;
      zi_heatflux[i] = 0.0;
      zi_lum_local[i] = 0.0;
   }  

  for(int i=0;i<nphi*ntheta*3;i++) {
    zi_heaterms[i] = 0.0;
    zi_heateave[i] = 0.0;
   }  

  return;
}
