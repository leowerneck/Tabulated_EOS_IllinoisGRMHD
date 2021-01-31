#ifndef __HARM_TABEOS_HELPERS__
#define __HARM_TABEOS_HELPERS__

// This is defined in the nuc_eos.hh file
#ifndef NTABLES
#define NTABLES (19)
#endif

//----------- EOS_Omni stuff -----------
// The following extern variables are declared in
// the EOS_Omni thorn file readtable.cc. They are
// initialized by the function nuc_eos_C_ReadTable(),
// in the same file.
namespace nuc_eos {
  extern double eos_rhomin,eos_yemin,eos_tempmin;
  extern double eos_rhomax,eos_yemax,eos_tempmax;
  extern double energy_shift;

  extern "C"
  void nuc_eos_m_kt1_press_eps_(const int *restrict n_in,
                                const double *restrict rho, 
                                const double *restrict temp,
                                const double *restrict ye,
                                double *restrict eps,
                                double *restrict prs,
                                int *restrict keyerr,
                                int *restrict anyerr);
}
namespace nuc_eos_private {
  extern int nrho,nye,ntemp;
  extern double drhoi,dtempi,dyei;
  extern double drhotempi,drhoyei,drhotempyei,dtempyei;
  extern double * restrict alltables;
  extern double * restrict logrho;
  extern double * restrict logtemp;
  extern double * restrict yes;
}

namespace HARM_TabEOS_helpers {

  static inline __attribute__((always_inline))
  int checkbounds_kt0_noTcheck(const double xrho, 
                               const double xye) {

    // keyerr codes:
    // 101 -- Y_e too high
    // 102 -- Y_e too low
    // 105 -- rho too high
    // 106 -- rho too low

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;

    if(xrho > eos_rhomax) {
      return 105;
    }
    if(xrho < eos_rhomin) {
      return 106;
    }
    if(xye > eos_yemax) {
      return 101;
    }
    if(xye < eos_yemin) {
      return 102;
    }
    return 0;
  }

  static inline __attribute__((always_inline))  
  void get_interp_spots(const double x,
                        const double y,
                        const double z,
                        double* restrict delx,
                        double* restrict dely,
                        double* restrict delz,
                        int* restrict idx) 
  {

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;

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

  static inline __attribute__((always_inline))  
  void get_interp_spots_d(const double x,
                          const double y,
                          const double z,
                          double* restrict delx,
                          double* restrict dely,
                          double* restrict delz,
                          int* restrict idx) 
  {

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;

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

  static inline __attribute__((always_inline))
  void nuc_eos_C_linterp_one(const int* restrict idx,
                             const double delx,
                             const double dely,
                             const double delz,
                             double* restrict f,
                             const int iv)
  {

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;

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

  static inline __attribute__((always_inline))
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

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;

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

  static inline __attribute__((always_inline))
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

  static inline __attribute__((always_inline))
  void bisection_enthalpy(const double lr, 
                          const double lt0,
                          const double ye,
                          double leps0,
                          const double prec,
                          double *restrict ltout,
                          const int iv,
                          int *restrict keyerrt) {

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;
    
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
    double lt1, lt2;
    double ltmin = logtemp[0];
    double ltmax = logtemp[ntemp-1];
    double f1,f2,ltmid,fmida_ent;
    double f1a[2] = {0};
    double f2a[2] = {0};
    double fmida[2] = {0};
    double delx,dely,delz;
    int idx[8];
 
    //  energy_shift = 0.0;

    // prepare
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

  static inline __attribute__((always_inline))
  void nuc_eos_findtemp_enthalpy(const double lr, 
                                 const double lt0,
                                 const double ye,
                                 const double enthalpyin,
                                 const double prec,
                                 double *restrict ltout,
                                 int *keyerrt) {

    // Leo says: added for easier interface with EOS_Omni
    using namespace nuc_eos;
    using namespace nuc_eos_private;

    // local variables
    const int itmax = 50; // use at most 10 iterations, then go to bisection
    double ldt;
    double ltn; // temp vars for temperature
    double oerr;
    const double ltmax = logtemp[ntemp-1]; // max temp
    const double ltmin = logtemp[0]; // min temp
    int it = 0;

    double dlenthalpydlt,denthalpydtc,dlpressdlt,dlepsdlt,depsdrho,dpdrho; // derivative dlogenthalpy/dlogT
    double lenthalpy,enthalpyc,lenthalpy0;
    double leps, lpress;
    // setting up some vars
    *keyerrt = 0;
    lenthalpy0 = enthalpyin;
    double lt = lt0;

    // step 1: do we already have the right temperature
    int idx[11];
    double delx,dely,delz;

    get_interp_spots_d(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,1,&leps,&depsdrho,&dlepsdlt,&dpdrho,&dlpressdlt);
    nuc_eos_C_linterp_one_d(idx,delx,dely,delz,0,&lpress,&depsdrho,&dlepsdlt,&dpdrho,&dlpressdlt);

    enthalpyc = exp(leps) - energy_shift + exp(lpress)/exp(lr);
    lenthalpy = enthalpyc;
    denthalpydtc = exp(leps)/exp(lt)*dlepsdlt + exp(lpress)/exp(lt)*dlpressdlt/exp(lr);
    dlenthalpydlt = exp(lt)/enthalpyc*denthalpydtc;

    if(fabs(lenthalpy-lenthalpy0) < prec*1.0e-2*fabs(lenthalpy0)) {
      *ltout = lt0;
      return;
    }

    oerr = 1.0e90;
    double fac = 1.0;

    while(it < itmax) {
      it++;

      // step 2: check if the two bounding values of the temperature
      //         give eps values that enclose the new eps.
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

      // drive the thing into the right direction
      double err = fabs(lenthalpy-lenthalpy0);
      if(oerr < err) fac *= 0.9;
      oerr = err;

      if(err < prec*1.0e-2*fabs(lenthalpy0)) {
	*ltout = lt;
	return;
      }
    } // while(it < itmax)

    bisection_enthalpy(lr,lt0,ye,lenthalpy0,prec,ltout,1,keyerrt);

    return;
  }

  
}

#endif // __HARM_TABEOS_HELPERS__
