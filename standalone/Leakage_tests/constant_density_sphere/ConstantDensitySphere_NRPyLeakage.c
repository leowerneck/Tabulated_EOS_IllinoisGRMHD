#include "Basic_defines.h"
#include "NRPyEOS.h"
#include "NRPyLeakage.h"

#ifdef IDX3D
#undef IDX3D
#endif
#define IDX3D(i0,i1,i2) ( (i0) + Nt0*( (i1) + Nt1*(i2) ) )

void dump_1d_data( const char *filename_prefix,
                   const int Nt0, const int Nt1, const int Nt2,
                   const int Ng0, const int Ng1, const int Ng2,
                   REAL **xx,
                   REAL **kappa_nue,
                   REAL **kappa_anue,
                   REAL **kappa_nux,
                   REAL **tau_nue,
                   REAL **tau_anue,
                   REAL **tau_nux ) {
  fprintf(stderr,"(ConstantDensitySphere) Dumping 1D data with prefix \"%s\"... ",filename_prefix);

  int midpoint[2];
  REAL x_midpoint[2];
  char filename[256];
  FILE *fp;
  // x-axis
  sprintf(filename,"%s_dump_x.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint  [0] = Nt1/2;
  midpoint  [1] = Nt2/2;
  x_midpoint[0] = xx[1][midpoint[0]];
  x_midpoint[1] = xx[2][midpoint[1]];
  for(int i0=Ng0;i0<Nt0-Ng0;i0++) {
    const REAL x0 = xx[0][i0];
    const int index = IDX3D(i0,midpoint[0],midpoint[1]);
    fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            x0,x_midpoint[0],x_midpoint[1],
            kappa_nue [0][index],kappa_nue [1][index],
            kappa_anue[0][index],kappa_anue[1][index],
            kappa_nux [0][index],kappa_nux [1][index],
            tau_nue   [0][index],tau_nue   [1][index],
            tau_anue  [0][index],tau_anue  [1][index],
            tau_nux   [0][index],tau_nux   [1][index]);
  }
  fclose(fp);

  // y-axis
  sprintf(filename,"%s_dump_y.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint  [0] = Nt0/2;
  midpoint  [1] = Nt2/2;
  x_midpoint[0] = xx[0][midpoint[0]];
  x_midpoint[1] = xx[2][midpoint[1]];
  for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
    const REAL x1 = xx[1][i1];
    const int index = IDX3D(midpoint[0],i1,midpoint[1]);
    fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            x_midpoint[0],x1,x_midpoint[1],
            kappa_nue [0][index],kappa_nue [1][index],
            kappa_anue[0][index],kappa_anue[1][index],
            kappa_nux [0][index],kappa_nux [1][index],
            tau_nue   [0][index],tau_nue   [1][index],
            tau_anue  [0][index],tau_anue  [1][index],
            tau_nux   [0][index],tau_nux   [1][index]);
  }
  fclose(fp);

  // z-axis
  sprintf(filename,"%s_dump_z.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint  [0] = Nt0/2;
  midpoint  [1] = Nt1/2;
  x_midpoint[0] = xx[0][midpoint[0]];
  x_midpoint[1] = xx[1][midpoint[1]];
  for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
    const REAL x2 = xx[2][i2];
    const int index = IDX3D(midpoint[0],midpoint[1],i2);
    fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
            x_midpoint[0],x_midpoint[1],x2,
            kappa_nue [0][index],kappa_nue [1][index],
            kappa_anue[0][index],kappa_anue[1][index],
            kappa_nux [0][index],kappa_nux [1][index],
            tau_nue   [0][index],tau_nue   [1][index],
            tau_anue  [0][index],tau_anue  [1][index],
            tau_nux   [0][index],tau_nux   [1][index]);
  }
  fclose(fp);
  fprintf(stderr,"done!\n");
}

void dump_2d_data( const char *filename_prefix,
                   const int Nt0, const int Nt1, const int Nt2,
                   const int Ng0, const int Ng1, const int Ng2,
                   REAL **xx,
                   REAL **kappa_nue,
                   REAL **kappa_anue,
                   REAL **kappa_nux,
                   REAL **tau_nue,
                   REAL **tau_anue,
                   REAL **tau_nux ) {

  fprintf(stderr,"(ConstantDensitySphere) Dumping 2D data with prefix \"%s\"... ",filename_prefix);

  int midpoint;
  REAL x_midpoint;
  char filename[256];
  FILE *fp;
  // xy-plane
  sprintf(filename,"%s_dump_xy.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint = Nt2/2;
  x_midpoint = xx[2][midpoint];
  for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
    const REAL x1 = xx[1][i1];
    for(int i0=Ng0;i0<Nt0-Ng0;i0++) {
      const REAL x0 = xx[0][i0];
      const int index = IDX3D(i0,i1,midpoint);
      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
              x0,x1,x_midpoint,
              kappa_nue [0][index],kappa_nue [1][index],
              kappa_anue[0][index],kappa_anue[1][index],
              kappa_nux [0][index],kappa_nux [1][index],
              tau_nue   [0][index],tau_nue   [1][index],
              tau_anue  [0][index],tau_anue  [1][index],
              tau_nux   [0][index],tau_nux   [1][index]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  // xz-plane
  sprintf(filename,"%s_dump_xz.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint = Nt1/2;
  x_midpoint = xx[1][midpoint];
  for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
    const REAL x2 = xx[2][i2];
    for(int i0=Ng0;i0<Nt0-Ng0;i0++) {
      const REAL x0 = xx[0][i0];
      const int index = IDX3D(i0,midpoint,i2);
      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
              x0,x_midpoint,x2,
              kappa_nue [0][index],kappa_nue [1][index],
              kappa_anue[0][index],kappa_anue[1][index],
              kappa_nux [0][index],kappa_nux [1][index],
              tau_nue   [0][index],tau_nue   [1][index],
              tau_anue  [0][index],tau_anue  [1][index],
              tau_nux   [0][index],tau_nux   [1][index]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  // yz-plane
  sprintf(filename,"%s_dump_yz.asc",filename_prefix);
  fp = fopen(filename,"w");
  midpoint = Nt0/2;
  x_midpoint = xx[0][midpoint];
  for(int i2=Ng2;i2<Nt2-Ng2;i2++) {
    const REAL x2 = xx[2][i2];
    for(int i1=Ng1;i1<Nt1-Ng1;i1++) {
      const REAL x1 = xx[1][i1];
      const int index = IDX3D(midpoint,i1,i2);
      fprintf(fp,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
              x_midpoint,x1,x2,
              kappa_nue [0][index],kappa_nue [1][index],
              kappa_anue[0][index],kappa_anue[1][index],
              kappa_nux [0][index],kappa_nux [1][index],
              tau_nue   [0][index],tau_nue   [1][index],
              tau_anue  [0][index],tau_anue  [1][index],
              tau_nux   [0][index],tau_nux   [1][index]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  fprintf(stderr,"done!\n");
}

void ConstantDensitySphere_NRPyLeakage(const NRPyEOS_params *restrict eos_params) {

  // Step 1: Set basic parameters
  const int N0     = 64;
  const int N1     = N0;
  const int N2     = N0;
  const int Ng0    =  2;
  const int Ng1    = Ng0;
  const int Ng2    = Ng0;
  const int Nt0    = N0 + 2*Ng0;
  const int Nt1    = N1 + 2*Ng1;
  const int Nt2    = N2 + 2*Ng2;
  const int Ntotal = Nt0 * Nt1 * Nt2;
  const REAL xmax  = +2;
  const REAL xmin  = -2;
  const REAL ymax  = xmax;
  const REAL ymin  = xmin;
  const REAL zmax  = xmax;
  const REAL zmin  = xmin;
  const REAL dx    = (xmax-xmin)/((REAL)N0);
  const REAL dy    = (xmax-xmin)/((REAL)N1);
  const REAL dz    = (xmax-xmin)/((REAL)N2);
  const REAL rSph  = 1;

  // Step 2: Allocate memory for the metric, opacities, and optical depths
  REAL **xx         = (REAL **)malloc(sizeof(REAL *)*3);
  xx[0]             = (REAL *)malloc(sizeof(REAL)*Nt0);
  xx[1]             = (REAL *)malloc(sizeof(REAL)*Nt1);
  xx[2]             = (REAL *)malloc(sizeof(REAL)*Nt2);
  REAL *gammaDD00   = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL *gammaDD11   = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL *gammaDD22   = (REAL *)malloc(sizeof(REAL)*Ntotal);
  REAL **kappa_nue  = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **kappa_anue = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **kappa_nux  = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_nue    = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_anue   = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_nux    = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_nue_p  = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_anue_p = (REAL **)malloc(sizeof(REAL *)*2);
  REAL **tau_nux_p  = (REAL **)malloc(sizeof(REAL *)*2);
  for(int i=0;i<2;i++) {
   kappa_nue [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   kappa_anue[i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   kappa_nux [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nue   [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_anue  [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nux   [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nue_p [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_anue_p[i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
   tau_nux_p [i] = (REAL *)malloc(sizeof(REAL)*Ntotal);
  }

  // Step 3: Define hydro quantities at the sphere interior and exterior
  //         The parameters below are extracted from https://arxiv.org/abs/2106.05356
  // Step 3.a: Interior of the sphere (sec. 4.2.1 of above referece)
  const REAL rho_interior = 9.8e13 * NRPyLeakage_units_cgs_to_geom_D;
  const REAL Y_e_interior = 0.1;
  const REAL T_interior   = 8.0;
  // Step 3.b: Exterior of the sphere (sec. 4.2.1 of above referece)
  const REAL rho_exterior = 6.0e7 * NRPyLeakage_units_cgs_to_geom_D;
  const REAL Y_e_exterior = 0.5;
  const REAL T_exterior   = 0.01;

  // Step 4: Compute opacities in the interior and exterior
  const REAL tau_nue_in [2] = {0.0,0.0};
  const REAL tau_anue_in[2] = {0.0,0.0};
  const REAL tau_nux_in [2] = {0.0,0.0};
  // Step 4.a: Interior
  REAL __attribute__((unused)) R_source,Q_source;
  REAL kappa_nue_interior[2], kappa_anue_interior[2], kappa_nux_interior[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(USE_NRPY_CONSTANTS,
                                                       eos_params,rho_interior,Y_e_interior,T_interior,
                                                       tau_nue_in,tau_anue_in,tau_nux_in,
                                                       &R_source,&Q_source,
                                                       kappa_nue_interior,kappa_anue_interior,kappa_nux_interior);
  // Step 4.b: Exterior
  REAL kappa_nue_exterior[2], kappa_anue_exterior[2], kappa_nux_exterior[2];
  NRPyLeakage_compute_GRMHD_source_terms_and_opacities(USE_NRPY_CONSTANTS,
                                                       eos_params,rho_exterior,Y_e_exterior,T_exterior,
                                                       tau_nue_in,tau_anue_in,tau_nux_in,
                                                       &R_source,&Q_source,
                                                       kappa_nue_exterior,kappa_anue_exterior,kappa_nux_exterior);

  // Step 5: Print basic information
  fprintf(stderr,"(ConstantDensitySphere) Test information:\n");
  fprintf(stderr,"(ConstantDensitySphere)     Domain properties:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - Sphere radius = %22.15e\n",rSph);
  fprintf(stderr,"(ConstantDensitySphere)         - xmin          = %22.15e\n",xmin);
  fprintf(stderr,"(ConstantDensitySphere)         - xmax          = %22.15e\n",xmax);
  fprintf(stderr,"(ConstantDensitySphere)         - ymin          = %22.15e\n",ymin);
  fprintf(stderr,"(ConstantDensitySphere)         - ymax          = %22.15e\n",ymax);
  fprintf(stderr,"(ConstantDensitySphere)         - zmin          = %22.15e\n",zmin);
  fprintf(stderr,"(ConstantDensitySphere)         - zmax          = %22.15e\n",zmax);
  fprintf(stderr,"(ConstantDensitySphere)         - Nx            = (%d) + 2x(%d)\n",N0,Ng0);
  fprintf(stderr,"(ConstantDensitySphere)         - Ny            = (%d) + 2x(%d)\n",N1,Ng1);
  fprintf(stderr,"(ConstantDensitySphere)         - Nz            = (%d) + 2x(%d)\n",N2,Ng2);
  fprintf(stderr,"(ConstantDensitySphere)         - dx            = %22.15e\n",dx);
  fprintf(stderr,"(ConstantDensitySphere)         - dy            = %22.15e\n",dy);
  fprintf(stderr,"(ConstantDensitySphere)         - dz            = %22.15e\n",dz);
  fprintf(stderr,"(ConstantDensitySphere)     Hydro quantities at sphere interior:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - rho           = %22.15e\n",rho_interior);
  fprintf(stderr,"(ConstantDensitySphere)         - Y_e           = %22.15e\n",Y_e_interior);
  fprintf(stderr,"(ConstantDensitySphere)         -  T            = %22.15e\n",T_interior);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nue   = %22.15e\n",kappa_nue_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nue   = %22.15e\n",kappa_nue_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_anue  = %22.15e\n",kappa_anue_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_anue  = %22.15e\n",kappa_anue_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nux   = %22.15e\n",kappa_nux_interior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nux   = %22.15e\n",kappa_nux_interior[1]);
  fprintf(stderr,"(ConstantDensitySphere)     Hydro quantities at sphere exterior:\n");
  fprintf(stderr,"(ConstantDensitySphere)         - rho           = %22.15e\n",rho_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         - Y_e           = %22.15e\n",Y_e_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         -  T            = %22.15e\n",T_exterior);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nue   = %22.15e\n",kappa_nue_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nue   = %22.15e\n",kappa_nue_exterior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_anue  = %22.15e\n",kappa_anue_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_anue  = %22.15e\n",kappa_anue_exterior[1]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_0_nux   = %22.15e\n",kappa_nux_exterior[0]);
  fprintf(stderr,"(ConstantDensitySphere)         - kappa_1_nux   = %22.15e\n",kappa_nux_exterior[1]);

  // Step 3: Loop over the grid, set opacities
  //         and initialize optical depth to zero
#pragma omp parallel for
  for(int i2=0;i2<Nt2;i2++) {
    const REAL z = zmin + (i2-Ng2+0.5)*dz;
    xx[2][i2] = z;
    for(int i1=0;i1<Nt1;i1++) {
      const REAL y = ymin + (i1-Ng1+0.5)*dy;
      xx[1][i1] = y;
      for(int i0=0;i0<Nt0;i0++) {
        const REAL x = xmin + (i0-Ng0+0.5)*dx;
        xx[0][i0] = x;

        // Step 3.a: Set local index
        const int index = IDX3D(i0,i1,i2);

        // Step 3.b: Initialize metric to flat space
        gammaDD00[index] = 1.0;
        gammaDD11[index] = 1.0;
        gammaDD22[index] = 1.0;

        // Step 3.c: Set local opacities
        const REAL r = sqrt( x*x + y*y + z*z );
        if( r > rSph ) {
          // Exterior
          kappa_nue [0][index] = kappa_nue_exterior [0];
          kappa_nue [1][index] = kappa_nue_exterior [1];
          kappa_anue[0][index] = kappa_anue_exterior[0];
          kappa_anue[1][index] = kappa_anue_exterior[1];
          kappa_nux [0][index] = kappa_nux_exterior [0];
          kappa_nux [1][index] = kappa_nux_exterior [1];
        }
        else {
          // Interior
          kappa_nue [0][index] = kappa_nue_interior [0];
          kappa_nue [1][index] = kappa_nue_interior [1];
          kappa_anue[0][index] = kappa_anue_interior[0];
          kappa_anue[1][index] = kappa_anue_interior[1];
          kappa_nux [0][index] = kappa_nux_interior [0];
          kappa_nux [1][index] = kappa_nux_interior [1];
        }

        // Step 3.d: Initialize optical depths to zero
        tau_nue   [0][index] = 0.0;
        tau_nue   [1][index] = 0.0;
        tau_anue  [0][index] = 0.0;
        tau_anue  [1][index] = 0.0;
        tau_nux   [0][index] = 0.0;
        tau_nux   [1][index] = 0.0;
      }
    }
  }

  // Step 4: Initial data dumps
  dump_1d_data("initial",Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,xx,kappa_nue,kappa_anue,kappa_nux,tau_nue,tau_anue,tau_nux);
  dump_2d_data("initial",Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,xx,kappa_nue,kappa_anue,kappa_nux,tau_nue,tau_anue,tau_nux);

  // Step 5: Now update the optical depth
  for(int n=0;n<20*N0;n++) {

    // Step 5.a: Copy optical depth to previous level
#pragma omp parallel for
    for(int i2=0;i2<Nt2;i2++) {
      for(int i1=0;i1<Nt1;i1++) {
        for(int i0=0;i0<Nt0;i0++) {
          const int index = IDX3D(i0,i1,i2);
          tau_nue_p [0][index] = tau_nue [0][index];
          tau_nue_p [1][index] = tau_nue [1][index];
          tau_anue_p[0][index] = tau_anue[0][index];
          tau_anue_p[1][index] = tau_anue[1][index];
          tau_nux_p [0][index] = tau_nux [0][index];
          tau_nux_p [1][index] = tau_nux [1][index];
        }
      }
    }

    // Step 5.b: Update optical depth
    NRPyLeakage_compute_optical_depths(Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,
                                       dx,dy,dz,
                                       gammaDD00,gammaDD11,gammaDD22,
                                       kappa_nue [0],kappa_nue [1],
                                       kappa_anue[0],kappa_anue[1],
                                       kappa_nux [0],kappa_nux [1],
                                       tau_nue_p [0],tau_nue_p [1],
                                       tau_anue_p[0],tau_anue_p[1],
                                       tau_nux_p [0],tau_nux_p [1],
                                       tau_nue   [0],tau_nue   [1],
                                       tau_anue  [0],tau_anue  [1],
                                       tau_nux   [0],tau_nux   [1]);
  }

  // Step 6: Final data dumps
  dump_1d_data("final",Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,xx,kappa_nue,kappa_anue,kappa_nux,tau_nue,tau_anue,tau_nux);
  dump_2d_data("final",Nt0,Nt1,Nt2,Ng0,Ng1,Ng2,xx,kappa_nue,kappa_anue,kappa_nux,tau_nue,tau_anue,tau_nux);

  // Step 7: Free memory
  for(int i=0;i<3;i++) free(xx[i]);
  free(xx);
  free(gammaDD00);
  free(gammaDD11);
  free(gammaDD22);
  for(int i=0;i<2;i++) {
    free(kappa_nue [i]);
    free(kappa_anue[i]);
    free(kappa_nux [i]);
    free(tau_nue   [i]);
    free(tau_anue  [i]);
    free(tau_nux   [i]);
    free(tau_nue_p [i]);
    free(tau_anue_p[i]);
    free(tau_nux_p [i]);
  }
  free(kappa_nue);
  free(kappa_anue);
  free(kappa_nux);
  free(tau_nue);
  free(tau_anue);
  free(tau_nux);
  free(tau_nue_p);
  free(tau_anue_p);
  free(tau_nux_p);
}
