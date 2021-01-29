/***********************************************************************************

Helper functions for the interpolation of EOS tables.
The functions are based on  Siegel, Mösta, Desai and Wu 2018, and Schneider, Roberts, Ott, 2017.

***********************************************************************************/

#include <math.h>
#include "nuc_eos.h"
#define DEBUG 0 
int checkbounds(const double xrho, 
		const double xtemp, 
		const double xye) {

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_rhomax) {

    fprintf(stderr,"rho: %15.6E eos_rhomax: %15.6E\n eos_rhomin: %15.6E\n",xrho,eos_rhomax,eos_rhomin);
    return 105;
  }
  if(xrho < eos_rhomin) {
    fprintf(stderr,"rho: %15.6E eos_rhomax: %15.6E\n eos_rhomin: %15.6E\n",xrho,eos_rhomax,eos_rhomin);

    return 106;
  }
  //printf("Inside helper rho");
  if(xye > eos_yemax) {
    fprintf(stderr,"ye: %15.6E eos_yemax: %15.6E\n ye_rhomin: %15.6E\n",xye,eos_yemax,eos_yemin);
    return 101;
  }
  if(xye < eos_yemin) {
    // this is probably not pure and should be removed
    fprintf(stderr,"xye: %15.6E eos_yemin: %15.6E\n",xye,eos_yemin);
    return 102;
  }
  //printf("Inside helper ye");
  //fprintf(stderr,"t: %15.6E eos_tmin: %15.6E\n eos_tmax: %15.6E\n",xtemp,eos_tempmax,eos_tempmin);
  if(xtemp > eos_tempmax) {
    fprintf(stderr,"t: %15.6E eos_tmin: %15.6E\n eos_tmax: %15.6E\n",xtemp,eos_tempmax,eos_tempmin);
    return 103;
  }
  if(xtemp < eos_tempmin) {
    fprintf(stderr, "temp_min %e\n", eos_tempmin);
    fprintf(stderr, "temp %e\n", xtemp);
    return 104;
  }
  //printf("Inside helper temp");
  return 0;
}

int checkbounds_kt0_noTcheck(const double xrho, 
		const double xye) {

	// keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 105 -- rho too high
  // 106 -- rho too low

//  printf("Inside cb, rhomi = %g\n",eos_rhomin);

  if(xrho > eos_rhomax) {
    fprintf(stderr,"rho: %15.6E eos_rhomax: %15.6E\n eos_rhomin: %15.6E\n",xrho,eos_rhomax,eos_rhomin);
    return 105;
  }
  if(xrho < eos_rhomin) {
    fprintf(stderr,"rho: %15.6E eos_rhomax: %15.6E\n eos_rhomin: %15.6E\n",xrho,eos_rhomax,eos_rhomin);
    return 106;
  }
//  printf("Inside cb, after rho\n");
  if(xye > eos_yemax) {
    fprintf(stderr,"ye: %15.6E eos_yemax: %15.6E\n eos_yemin: %15.6E\n",xye,eos_yemax,eos_yemin);
    return 101;
  }
  if(xye < eos_yemin) {
    return 102;
  }
//  printf("Inside cb, after ye\n");
  return 0;
}

void get_interp_spots(const double x,
		      const double y,
		      const double z,
		      double* restrict delx,
		      double* restrict dely,
		      double* restrict delz,
		      int* restrict idx) 
{

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1 + (int)( (y - logtemp[0] - 1.0e-10) * dtempi );
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iy = MAX( 1, MIN( iy, ntemp-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = logtemp[iy] - y;
  *delz = yes[iz] - z;

  return;
}


void get_interp_spots_d(const double x,
		      const double y,
		      const double z,
		      double* restrict delx,
		      double* restrict dely,
		      double* restrict delz,
		      int* restrict idx) 
{

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1 + (int)( (y - logtemp[0] - 1.0e-10) * dtempi );
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iy = MAX( 1, MIN( iy, ntemp-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));
  idx[8] = NTABLES*((ix+1) + nrho*(iy + ntemp*iz));
  idx[9] = NTABLES*(ix + nrho*((iy+1) + ntemp*iz));
  idx[10] = NTABLES*(ix + nrho*(iy + ntemp*(iz+1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = logtemp[iy] - y;
  *delz = yes[iz] - z;

  return;
}

void get_interp_spots_linT_low(const double x,
			       const double y,
			       const double z,
			       double* restrict delx,
			       double* restrict dely,
			       double* restrict delz,
			       int* restrict idx) 
{

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}

void get_interp_spots_linT_low_eps(const double x,
				   const double y,
				   const double z,
				   double* restrict delx,
				   double* restrict dely,
				   double* restrict delz,
				   int* restrict idx) 
{

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = (ix + nrho*(iy + ntemp*iz));
  idx[1] = ((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = (ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = (ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = ((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = ((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = (ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = ((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}




void get_interp_spots_linT_low_entropy(const double x,
           const double y,
           const double z,
           double* restrict delx,
           double* restrict dely,
           double* restrict delz,
           int* restrict idx) 
{

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = (ix + nrho*(iy + ntemp*iz));
  idx[1] = ((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = (ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = (ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = ((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = ((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = (ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = ((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}

void nuc_eos_C_linterp_one_d(const int* restrict idx, 
			   const double delx, 
			   const double dely,
			   const double delz,
			   const int iv,
			   double* restrict f,
                           double* restrict depsdrho, 
			   double* restrict depsdt, 
  			   double* restrict dPdrho, 
			   double* restrict dPdt)
{

  // helper variables
  double fh[11], a[11];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];
  fh[8] = alltables[iv+idx[8]];
  fh[9] = alltables[iv+idx[9]];
  fh[10] = alltables[iv+idx[10]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dtempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drhotempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dtempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );
  a[8] = 0.5 * drhoi *   ( fh[1] - fh[8] );
  a[9] = 0.5 * dtempi *   ( fh[2] - fh[9] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  if (iv==0) {
    *dPdt = -a[2];
    *dPdrho = -a[1];
  }
  if (iv==1) {
    *depsdt = -a[2];
    *depsdrho = -a[1];
  }
  return;
}


void nuc_eos_C_linterp_one(const int* restrict idx, 
			   const double delx, 
			   const double dely,
			   const double delz,
			   double* restrict f,
			   const int iv)
{

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dtempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drhotempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dtempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}

void nuc_eos_C_linterp_one_linT_low(const int* restrict idx, 
				    const double delx, 
				    const double dely,
				    const double delz,
				    double* restrict f,
				    const int iv)
{

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}


void nuc_eos_C_linterp_one_linT_low_eps(const int* restrict idx, 
					const double delx, 
					const double dely,
					const double delz,
					double* restrict f)

{

  // helper variables
  double fh[8], a[8];

  fh[0] = epstable[idx[0]];
  fh[1] = epstable[idx[1]];
  fh[2] = epstable[idx[2]];
  fh[3] = epstable[idx[3]];
  fh[4] = epstable[idx[4]];
  fh[5] = epstable[idx[5]];
  fh[6] = epstable[idx[6]];
  fh[7] = epstable[idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}




double linterp2D(const double *restrict xs, 
		 const double *restrict ys, 
		 const double *restrict fs, 
		 const double x, 
		 const double y)
{

  //  2     3 
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid
  
  double dxi = 1./(xs[1]-xs[0]);
  double dyi = 1./(ys[1]-ys[0]); // x*1./y uses faster instructions than x/y
  double t1 = (fs[1]-fs[0])*dxi * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])*dxi * (x - xs[0]) + fs[2];

  return (t2 - t1)*dyi * (y-ys[0]) + t1;
}

void bisection_enthalpy(const double lr, 
	       const double lt0,
	       const double ye,
	       double leps0,
	       const double prec,
	       double *restrict ltout,
	       const int iv,
	       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on


  int bcount = 0; 
  int bcount2 = 0; 
  int maxbcount = 160;
  int itmax = 200;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp = log(1.2);
  const double dltm = log(0.8);

  const double leps0_prec = leps0*prec;

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid,fmida_ent;
  double f1a[2] = {0};
  double f2a[2] = {0};
  double fmida[2] = {0};
  double delx,dely,delz;
  int idx[8];
 
//  energy_shift = 0.0;

  // prepare
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);  
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);  
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);


  double f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
  double f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
  f1=f1a_ent-leps0;
  f2=f2a_ent-leps0;
#if DEBUG
  printf("Bisection start: Count = %d, press = %g, eps = %g, f1a_ent = %g, f2a_ent = %g, enthalpy_in = %g\n", bcount,f1a[0], f1a[1], f1a_ent, f2a_ent, leps0);
  printf("Count = %d, f1a_ent = %g, f2a_ent = %g, f1 = %g, f2 = %g, ltmin=%g\n", bcount,f1a_ent,f2a_ent,f1,f2,ltmin);
#endif
  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);

    get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);  
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);  
    
    get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0);  
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1);

    f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
    f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
    f1=f1a_ent-leps0;
    f2=f2a_ent-leps0;

#if DEBUG
    printf("Count = %d, press = %g, eps = %g, enthalpy_in = %g\n", bcount,f1a[0], f1a[1], leps0);
    printf("Count = %d, press2 = %g, eps2 = %g, enthalpy_in = %g\n", bcount,f2a[0], f2a[1], leps0);
    printf("Count = %d, lt1 = %g, lt2 = %g\n", bcount,lt1,lt2);
    printf("Count = %d, f1a_ent = %g, f2a_ent = %g, f1 = %g, f2 = %g\n", bcount,f1a_ent,f2a_ent,f1,f2);
#endif

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a_ent,f2a_ent,leps0);
#endif
    bcount++;
    if(bcount >= maxbcount && bcount2 == 0) {
      bcount2 = 1;
      bcount = 0;
      get_interp_spots(lr,ltmin,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0);  
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1);

      double ent_min = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);
      leps0 = fmax(ent_min*0.99,ent_min*1.01);
#if DEBUG
      printf("Adjusting ent_min: Old ent = %g\n", leps0);
      printf("Adjusting ent_min: New ent = %g\n", leps0);
#endif      
    }
    else if(bcount >= maxbcount && bcount2 == 1) { 
      *keyerrt = 667;
      return;
    }
    else {
    }

  } // while

#if DEBUG
    fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    f2,lt,dlt);
    fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[0],0); 
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a[1],1); 
  f1a_ent = exp(f1a[1])-energy_shift+exp(f1a[0])/exp(lr);

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[0],0); 
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a[1],1); 

  f2a_ent = exp(f2a[1])-energy_shift+exp(f2a[0])/exp(lr);
  int it;
  for(it=0;it<itmax;it++) {
    ltmid = log((exp(lt1) + exp(lt2))*0.5);

    get_interp_spots(lr,ltmid,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&fmida[0],0); 
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&fmida[1],1);
    fmida_ent = exp(fmida[1])-energy_shift+exp(fmida[0])/exp(lr);

    if ((f1a_ent - leps0)*(fmida_ent-leps0) < 0.0) {
	    f2a_ent = fmida_ent;
	    lt2 = ltmid;
    } 
    else if ((f2a_ent - leps0)*(fmida_ent-leps0) < 0.0) {
	    f1a_ent = fmida_ent;
	    lt1 = ltmid;
    } 
    else {
#if DEBUG
      printf("Bracketing failed for bisection interval\n");
#endif
    }  

#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E, crit = %15.6E,leps0_prec = %15.6E\n",
	    it,fmid,f2a_ent,lt,ltmid,dlt,fabs(leps0-f2a_ent),leps0_prec);
#endif

#if DEBUG
     printf("It = %d, press2 = %g, eps2 = %g, enthalpy = %g, enthalpy_in = %g\n", it,f2a[0], f2a[1], f2a_ent, leps0);
     printf("It = %d, ltmid = %g, lt0 = %g\n", it,ltmid,lt0);
     printf("It = %d, f2a_ent = %g, err = %g\n", it,f2a_ent,fabs(1.0-f2a_ent/leps0));
#endif

    if(fabs(leps0-fmida_ent) <= fabs(leps0_prec*5.0e-2)) {

#if DEBUG
      printf("Bisection found correct temp: It = %d, f2a_ent = %g, err = %g, threshold = %g\n", it,f2a_ent,fabs(leps0-f2a_ent),fabs(leps0_prec*1.0e-2));
#endif
      *ltout = ltmid;
      return;
    } 
  } // for it = 0
  *keyerrt = 667;
  return;
} // bisection


void bisection(const double lr, 
	       const double lt0,
	       const double ye,
	       double leps0,
	       const double prec,
	       double *restrict ltout,
	       const int iv,
	       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on


  int bcount = 0; 
  int maxbcount = 80;
  int itmax = 100;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp = log(1.2);
  const double dltm = log(0.8);

  const double leps0_prec = leps0*prec;

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;
  double delx,dely,delz;
  int idx[8];
#if DEBUG  
  printf("ltmin = %g\n", ltmin); 
#endif
  // prepare
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);


  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);

  f1=f1a-leps0;
  f2=f2a-leps0;

  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);
    get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  
    
    get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);

    f1=f1a-leps0;
    f2=f2a-leps0;

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a,f2a,leps0);
#endif

    bcount++;
    if(bcount >= maxbcount) {
      *keyerrt = 667;
      return;
    }
  } // while

  if(f1 < 0.0) {
    lt = lt1;
    dlt = lt2 - lt1;
  } else {
    lt = lt2;
    dlt = lt1 - lt2;
  }

#if DEBUG
    fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    f2,lt,dlt);
    fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  int it;
  for(it=0;it<itmax;it++) {
    dlt = dlt * 0.5;
    ltmid = lt + dlt;
    get_interp_spots(lr,ltmid,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);
    
    fmid=f2a-leps0;
    if(fmid <= 0.0) lt=ltmid;
#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E, crit = %15.6E,leps0_prec = %15.6E\n",
	    it,fmid,f2a,lt,ltmid,dlt,fabs(leps0-f2a),leps0_prec);
#endif

    if(fabs(leps0-f2a) <= fabs(leps0_prec)) {
      *ltout = ltmid;
      return;
    }
  } // for it = 0

  *keyerrt = 667;
  return;
} // bisection

void nuc_eos_findtemp_enthalpy(const double lr, 
		      const double lt0,
		      const double ye,
		      const double enthalpyin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

    // local variables
    const int itmax = 50; // use at most 10 iterations, then go to bisection
    double dentdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double ent,ent0; // temp vars for eps
    double ltn, lt, lt1; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    const int iv = 1;
    int it = 0;

    double dlenthalpydlt,denthalpydtc,dlenthalpydltc,dlpressdlt,dlepsdlt,depsdrho,dpdrho; // derivative dlogenthalpy/dlogT
    double lenthalpy,enthalpy,lenthalpyc,enthalpyc,lenthalpy0,enthalpy0,lenthalpy1,enthalpy1;
    double leps, lpress;
    double presst1, presst2, epst1, epst2, enthalpyt1,enthalpyt2,enthalpyc1,enthalpyc2;
    // setting up some vars
    *keyerrt = 0;
    lenthalpy0 = enthalpyin;
    lenthalpy1 = lenthalpy0;
    lt = lt0;
    lt1 = lt;

    // step 1: do we already have the right temperature
    int idx[11];
    double delx,dely,delz;

    get_interp_spots_d(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,1,&leps,&depsdrho,&dlepsdlt,&dpdrho,&dlepsdlt);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,0,&lpress,&depsdrho,&dlepsdlt,&dpdrho,&dlpressdlt);


    enthalpyc = exp(leps) - energy_shift + exp(lpress)/exp(lr);
    lenthalpy = enthalpyc;
    denthalpydtc = exp(leps)/exp(lt)*dlepsdlt + exp(lpress)/exp(lt)*dlpressdlt/exp(lr);
    dlenthalpydlt = exp(lt)/enthalpyc*denthalpydtc;



    if(fabs(lenthalpy-lenthalpy0) < prec*1.0e-2*fabs(lenthalpy0)) {
      *ltout = lt0;
      return;
    }

    lt1 = lt;
    lenthalpy1 = lenthalpy;


    //fprintf(stderr,"it: %d t: %15.6E lenthalpy0: %15.6E lenthalpy: %15.6E del: %15.6E\n",
    //    it,lt, lenthalpy0,lenthalpy,fabs(ent-ent0)/(fabs(ent0)));
    
    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

      double enthalpy1, enthalpy2;
      // lower temperature
      {
	// get data at 4 points
	double fs[4],fs1[4];
	// point 0
	int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	int ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	fs1[0] = alltables[ifs1];
	// point 1 
	ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	fs1[1] = alltables[ifs1];
	// point 2 
	ifs = 0 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	fs1[2] = alltables[ifs1];
	// point 3
	ifs = 0 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	fs[3] = alltables[ifs];
	fs1[3] = alltables[ifs1];

	presst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
 	epst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs1, lr, ye);
        enthalpyc1 =  exp(epst1) - energy_shift + exp(presst1)/exp(lr);


      }
      // upper temperature
      {
	// get data at 4 points
	double fs[4],fs1[4];
	// point 0
	int ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	int ifs1 = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	fs1[0] = alltables[ifs1];
	// point 1 
	ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	fs1[1] = alltables[ifs1];
	// point 2 
	ifs = 0 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	fs1[2] = alltables[ifs1];
	// point 3
	ifs = 0 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	ifs1 = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	fs[3] = alltables[ifs];
	fs1[3] = alltables[ifs1];
      
 	presst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
 	epst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs1, lr, ye);
        enthalpyc2 = exp(epst2) - energy_shift + exp(presst2)/exp(lr);
      }

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
//      if((lenthalpy0 - enthalpyc1) * (lenthalpy0 - enthalpyc2) <= 0.) {
//      
//	*ltout = (logtemp[itemp]-logtemp[itemp-1]) / (enthalpyc2 - enthalpyc1) * 
//	  (lenthalpy0 - enthalpyc1) + logtemp[itemp-1];
//#if DEBUG
//        printf("Found correct temp in bracketing\n", it);
//#endif
//	return;
//      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      denthalpydtc = exp(leps)/exp(lt)*dlepsdlt + exp(lpress)/exp(lt)*dlpressdlt/exp(lr);
      dlenthalpydlt = exp(lt)*denthalpydtc;
      
      ldt = -(lenthalpy - lenthalpy0) / dlenthalpydlt;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,1);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&lpress,0);
      enthalpyc = exp(leps) - energy_shift + exp(lpress)/exp(lr);
      lenthalpy = enthalpyc;
      //fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      //fprintf(stderr,"it: %d t: %15.6E lenthalpy0: %15.6E lenthalpy: %15.6E del: %15.6E\n",
       // it,lt, lenthalpy0,lenthalpy,fabs(ent-ent0)/(fabs(ent0)));

#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
	      it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
      // drive the thing into the right direction
      double err = fabs(lenthalpy-lenthalpy0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*1.0e-2*fabs(lenthalpy0)) {
	*ltout = lt;
#if DEBUG
        printf("Found correct temp in %d iterations\n", it);
#endif
	return;
      }
    } // while(it < itmax)

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection_enthalpy(lr,lt0,ye,lenthalpy0,prec,ltout,1,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
    }      
#endif


    return;
}


void nuc_eos_findtemp(const double lr, 
		      const double lt0,
		      const double ye,
		      const double entin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

    // local variables
    const int itmax = 30; // use at most 10 iterations, then go to bisection
    double dentdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double ent,ent0; // temp vars for eps
    double ltn, lt; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    const int iv = 1;
    int it = 0;

    // setting up some vars
    *keyerrt = 0;
    ent0 = entin;
    lt = lt0;


    // step 1: do we already have the right temperature
    int idx[8];
    double delx,dely,delz;
    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);

    //fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
    //it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
    //fprintf(stderr, "ltmin %e, ltmax %e\n", ltmin,ltmax);

#if DEBUG
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif

    if(fabs(ent-ent0) < prec*fabs(ent0)) {
      *ltout = lt0;       
      return;
    }

    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

      double ent1, ent2;
      // lower temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	ent1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }
      // upper temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	ent2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }

      //fprintf(stderr, "ent0 %e, ent1 %e, ent2 %e, ent %e\n", ent0,ent1,ent2,ent);

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
      if((ent0 - ent1) * (ent0 - ent2) <= 0.) {
      
	*ltout = (logtemp[itemp]-logtemp[itemp-1]) / (ent2 - ent1) * 
	  (ent0 - ent1) + logtemp[itemp-1];    
	return;
      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      dentdlti = (logtemp[itemp]-logtemp[itemp-1])/(ent2-ent1);
      ldt = -(ent - ent0) * dentdlti * fac;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);
      
      //fprintf(stderr,"irho %d itemp %d iye %d\n",irho,itemp,iye);
      //fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
       // it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
      //fprintf(stderr, "ltim %e ltmax %e \n", ltmin,ltmax);
#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
	      it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
      // drive the thing into the right direction
      double err = fabs(ent-ent0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*fabs(ent0)) {
	*ltout = lt;     
	return;
      }

    } // while(it < itmax)

#if DEBUG
#endif

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    //fprintf(stderr, "lr %e,lt0 %e, ent0 %e, ltout %e, keyerrt %i\n", lr,lt0, ent0, ltout, *keyerrt);
    bisection(lr,lt0,ye,ent0,prec,ltout,iv,keyerrt);
//#if DEBUG
    if(*keyerrt==667) {
     // fprintf(stderr,"This is worse. Bisection failed!\n");
    }      
//#endif


    return;
  }



void nuc_eos_findtemp_entropy(const double lr, 
          const double lt0,
          const double ye,
          const double entin,
          const double prec,
          double *restrict ltout,
          int *keyerrt) {

    // local variables
    const int itmax = 30; // use at most 10 iterations, then go to bisection
    double dentdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double ent,ent0; // temp vars for entropy
    double ltn, lt; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    const int iv = 2;
    int it = 0;

    // setting up some vars
    *keyerrt = 0;
    ent0 = entin;
    lt = lt0;


    // step 1: do we already have the right temperature
    int idx[8];
    double delx,dely,delz;
    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);

    //fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
    //it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
    //fprintf(stderr, "ltmin %e, ltmax %e\n", ltmin,ltmax);

#if DEBUG
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
      it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif

    if(fabs(ent-ent0) < prec*fabs(ent0)) {
      *ltout = lt0;       
      return;
    }

    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give ent values that enclose the new ent.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

      double ent1, ent2;
      // lower temperature
      {
  // get data at 4 points
  double fs[4];
  // point 0
  int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
  fs[0] = alltables[ifs];
  // point 1 
  ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
  fs[1] = alltables[ifs];
  // point 2 
  ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
  fs[2] = alltables[ifs];
  // point 3
  ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
  fs[3] = alltables[ifs];
      
  ent1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }
      // upper temperature
      {
  // get data at 4 points
  double fs[4];
  // point 0
  int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
  fs[0] = alltables[ifs];
  // point 1 
  ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
  fs[1] = alltables[ifs];
  // point 2 
  ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
  fs[2] = alltables[ifs];
  // point 3
  ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
  fs[3] = alltables[ifs];
      
  ent2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }

      //fprintf(stderr, "ent0 %e, ent1 %e, ent2 %e, ent %e\n", ent0,ent1,ent2,ent);

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
      if((ent0 - ent1) * (ent0 - ent2) <= 0.) {
      
  *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (ent2 - ent1) * 
    (ent0 - ent1) + logtemp[itemp-1];    
  return;
      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      dentdlti = (logtemp[itemp]-logtemp[itemp-1])/(ent2-ent1);
      ldt = -(ent - ent0) * dentdlti * fac;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);
      
      //fprintf(stderr,"irho %d itemp %d iye %d\n",irho,itemp,iye);
      //fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
       // it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
      //fprintf(stderr, "ltim %e ltmax %e \n", ltmin,ltmax);
#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
        it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
      // drive the thing into the right direction
      double err = fabs(ent-ent0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*fabs(ent0)) {
  *ltout = lt;     
  return;
      }

    } // while(it < itmax)

#if DEBUG
#endif

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    //fprintf(stderr, "lr %e,lt0 %e, ent0 %e, ltout %e, keyerrt %i\n", lr,lt0, ent0, ltout, *keyerrt);
    bisection(lr,lt0,ye,ent0,prec,ltout,iv,keyerrt);
//#if DEBUG
    if(*keyerrt==667) {
     // fprintf(stderr,"This is worse. Bisection failed!\n");
    }      
//#endif


    return;
  }



void nuc_eos_findtemp_press(const double lr, 
				  const double lt0,
				  const double ye,
				  const double pressin,
				  const double prec,
				  double *restrict ltout,
				  int *keyerrt) {

    // local variables
    const int itmax = 200; // use at most 10 iterations, then go to bisection
    double dpressdlti; // 1 / derivative dentropy/dlogT
    double ldt;
    double press,press0; // temp vars for eps
    double ltn, lt; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    const int iv = 2;
    int it = 0;

    // setting up some vars
    *keyerrt = 0;
    press0 = pressin;
    lt = lt0;

    // step 1: do we already have the right temperature
    int idx[8];
    double delx,dely,delz;
    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&press,iv);
#if DEBUG
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt,press,press0,fabs(press-press0)/(fabs(press0)));
#endif
    if(fabs(press-press0) < prec*fabs(press0)) {
      *ltout = lt0;
      return;
    }

    oerr = 1.0e90;
    double fac = 1.0;
    const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
    const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
      const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

      double press1, press2;
      // lower temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	press1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }
      // upper temperature
      {
	// get data at 4 points
	double fs[4];
	// point 0
	int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
	fs[0] = alltables[ifs];
	// point 1 
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
	fs[1] = alltables[ifs];
	// point 2 
	ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
	fs[2] = alltables[ifs];
	// point 3
	ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
	fs[3] = alltables[ifs];
      
	press2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
      }

      // Check if we are already bracketing the input internal
      // energy. If so, interpolate for new T.
      if((press0 - press1) * (press0 - press2) <= 0.) {
      
	*ltout = (logtemp[itemp]-logtemp[itemp-1]) / (press2 - press1) * 
	  (press0 - press1) + logtemp[itemp-1];
     
	return;
      }

      // well, then do a Newton-Raphson step
      // first, guess the derivative
      dpressdlti = (logtemp[itemp]-logtemp[itemp-1])/(press2-press1);
      ldt = -(press - press0) * dpressdlti * fac;

      ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
      lt = ltn;

      get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&press,iv);

#if DEBUG
      fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
      fprintf(stderr,"it: %d t: %15.6E press: %15.6E press0: %15.6E del: %15.6E\n",
	      it,lt, press,press0,fabs(press-press0)/(fabs(press0)));
#endif
      // drive the thing into the right direction
      double err = fabs(press-press0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*fabs(press0)) {
	*ltout = lt;
	return;
      }



    } // while(it < itmax)

    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection(lr,lt0,ye,press0,prec,ltout,2,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
      abort();
    }      
#endif


    return;
  }
