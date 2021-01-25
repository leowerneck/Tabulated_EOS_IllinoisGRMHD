int con2prim( const igm_eos_parameters eos,
              const int index,const int i,const int j,const int k,
              const CCTK_REAL *restrict X,const CCTK_REAL *restrict Y,const CCTK_REAL *restrict Z,
              const CCTK_REAL *restrict METRIC,const CCTK_REAL *restrict METRIC_PHYS,const CCTK_REAL *restrict METRIC_LAP_PSI4,
              const CCTK_REAL g4dn[NDIM][NDIM],const CCTK_REAL g4up[NDIM][NDIM],
              CCTK_REAL *restrict CONSERVS,CCTK_REAL *restrict PRIMS,
              output_stats& stats ) {

  // declare some variables for HARM.
  CCTK_REAL cons[numcons];
  CCTK_REAL prim[numprims];

  /*
    -- Driver for new prim. var. solver.  The driver just translates
    between the two sets of definitions for U and P.  The user may
    wish to alter the translation as they see fit.


    //        /     rho u^t     \                           //
    //    U = | T^t_t + rho u^t | * sqrt(-det(g_{\mu\nu}))  //
    //        |     T^t_i       |                           //
    //        \      B^i        /                           //
    //                                                      //
    //        /     rho     \                               //
    //    P = |     uu      |                               //
    //        | \tilde{u}^i |                               //
    //        \     B^i     /                               //

    (above equations have been fixed by Yuk Tung & Zach)
  */

  // U[NPR]    = conserved variables (current values on input/output);
  // g4dn[NDIM][NDIM] = covariant form of the 4-metric ;
  // g4up[NDIM][NDIM] = contravariant form of the 4-metric ;
  // gdet             = sqrt( - determinant of the 4-metric) ;
  // prim[NPR] = primitive variables (guess on input, calculated values on
  //                     output if there are no problems);

  // U[1]   =
  // U[2-4] =  stildei + rhostar

  CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR  ];
  CCTK_REAL mhd_st_x_orig = CONSERVS[STILDEX  ];
  CCTK_REAL mhd_st_y_orig = CONSERVS[STILDEY  ];
  CCTK_REAL mhd_st_z_orig = CONSERVS[STILDEZ  ];
  CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];


  // Other ideas for setting the gamma speed limit
  //CCTK_REAL GAMMA_SPEED_LIMIT = 100.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=500.0;
  //if(METRIC_LAP_PSI4[PSI6]>Psi6threshold) GAMMA_SPEED_LIMIT=100.0;

  //FIXME: Only works if poisoning is turned on. Otherwise will access unknown memory. This trick alone speeds up the whole code (Cowling) by 2%.
  //int startguess=0;
  //if(std::isnan(PRIMS[VX])) startguess=1;
  int startguess=1;

  CCTK_REAL u0L=1.0;
  CCTK_REAL K_ppoly_tab,Gamma_ppoly_tab;

  for(int which_guess=startguess;which_guess<3;which_guess++) {

    // Set the conserved variables required by the con2prim routine
    set_cons_from_PRIMS_and_CONSERVS( eos, METRIC,METRIC_LAP_PSI4,PRIMS,CONSERVS, cons );

    // Set primitive guesses
    set_prim_from_PRIMS_and_CONSERVS( eos, which_guess,METRIC,METRIC_LAP_PSI4,PRIMS,CONSERVS, cons,prim );

    /************* Conservative-to-primitive recovery ************/
    int check = con2prim_Noble2D(eos,g4dn,g4up,cons,prim);
    /*************************************************************/

    // Use the new Font fix subroutine
    int font_fix_applied=0;
    if(check!=0) {
      font_fix_applied=1;
      CCTK_REAL u_xl=1e100, u_yl=1e100, u_zl=1e100; // Set to insane values to ensure they are overwritten.
      /************************
       * New Font fix routine *
       ************************/
      check = font_fix__hybrid_EOS(eos,METRIC_PHYS,METRIC_LAP_PSI4,CONSERVS,PRIMS, u_xl,u_yl,u_zl);

      //Translate to HARM primitive now:
      prim[UTCON1] = METRIC_PHYS[GUPXX]*u_xl + METRIC_PHYS[GUPXY]*u_yl + METRIC_PHYS[GUPXZ]*u_zl;
      prim[UTCON2] = METRIC_PHYS[GUPXY]*u_xl + METRIC_PHYS[GUPYY]*u_yl + METRIC_PHYS[GUPYZ]*u_zl;
      prim[UTCON3] = METRIC_PHYS[GUPXZ]*u_xl + METRIC_PHYS[GUPYZ]*u_yl + METRIC_PHYS[GUPZZ]*u_zl;
      if (check==1) {
        CCTK_VInfo(CCTK_THORNSTRING,"Font fix failed!");
        CCTK_VInfo(CCTK_THORNSTRING,"i,j,k = %d %d %d, stats.failure_checker = %d x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",i,j,k,stats.failure_checker,X[index],Y[index],Z[index],index,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
      }
    }
    stats.failure_checker+=font_fix_applied*10000;
    stats.font_fixed=font_fix_applied;
    /*************************************************************/


    if(check==0) {
      //Now that we have found some solution, we first limit velocity:
      //FIXME: Probably want to use exactly the same velocity limiter function here as in mhdflux.C
      CCTK_REAL utx_new = prim[UTCON1];
      CCTK_REAL uty_new = prim[UTCON2];
      CCTK_REAL utz_new = prim[UTCON3];

      //Velocity limiter:
      CCTK_REAL gijuiuj = METRIC_PHYS[GXX]*SQR(utx_new ) +
        2.0*METRIC_PHYS[GXY]*utx_new*uty_new + 2.0*METRIC_PHYS[GXZ]*utx_new*utz_new +
        METRIC_PHYS[GYY]*SQR(uty_new) + 2.0*METRIC_PHYS[GYZ]*uty_new*utz_new +
        METRIC_PHYS[GZZ]*SQR(utz_new);
      CCTK_REAL au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
      u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];

      // *** Limit velocity
      stats.vel_limited=0;
      if (au0m1 > 0.9999999*(eos.W_max-1.0)) {
        CCTK_REAL fac = sqrt((SQR(eos.W_max)-1.0)/(SQR(1.0+au0m1) - 1.0));
        utx_new *= fac;
        uty_new *= fac;
        utz_new *= fac;
        gijuiuj = gijuiuj * SQR(fac);
        au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
        // Reset rho_b and u0
        u0L = (au0m1+1.0)*METRIC_LAP_PSI4[LAPSEINV];
        prim[RHO] =  rho_star_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        stats.vel_limited=1;
        stats.failure_checker+=1000;
      } //Finished limiting velocity


      //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.
      if(font_fix_applied==1) {
        prim[RHO] = rho_star_orig/(METRIC_LAP_PSI4[LAPSE]*u0L*METRIC_LAP_PSI4[PSI6]);
        //Next set P = P_cold:
        CCTK_REAL P_cold;

        /**********************************
         * Piecewise Polytropic EOS Patch *
         *  Finding Gamma_ppoly_tab and K_ppoly_tab *
         **********************************/
        /* Here we use our newly implemented
         * find_polytropic_K_and_Gamma() function
         * to determine the relevant polytropic
         * Gamma and K parameters to be used
         * within this function.
         */
        int polytropic_index = find_polytropic_K_and_Gamma_index(eos,prim[RHO]);
        K_ppoly_tab     = eos.K_ppoly_tab[polytropic_index];
        Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];

        // After that, we compute P_cold
        P_cold = K_ppoly_tab*pow(prim[RHO],Gamma_ppoly_tab);

        prim[UU] = P_cold/(Gamma_ppoly_tab-1.0);
      } //Finished setting remaining primitives if there was a Font fix.

      /* Set rho_b */
      PRIMS[RHOB] = prim[RHO];

      /***************
       * PPEOS Patch *
       * Hybrid EOS  *
       ***************
       */
      /* We now compute the pressure as a function
       * of rhob, P_cold, eps_cold, and u = rhob*eps,
       * using the function pressure_rho0_u(), which
       * implements the equation:
       * .-------------------------------------------------------------.
       * | p(rho_b,u) = P_cold + (Gamma_th - 1)*(u - rho_b * eps_cold) |
       * .-------------------------------------------------------------.
       */
      PRIMS[PRESSURE] = pressure_rho0_u(eos, prim[RHO],prim[UU]);

      /* Already set u0L. */
      PRIMS[VX]       = utx_new/u0L - METRIC[SHIFTX];
      PRIMS[VY]       = uty_new/u0L - METRIC[SHIFTY];
      PRIMS[VZ]       = utz_new/u0L - METRIC[SHIFTZ];

      return 0;
      
    } else {
      //If we didn't find a root, then try again with a different guess.
    }
  }
  CCTK_VInfo(CCTK_THORNSTRING,"Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rho_b_atm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e",
	     tau_orig,rho_star_orig,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig/METRIC_LAP_PSI4[PSI6],eos.rho_atm,PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[LAPSE]);
  return 1;
}

#include "eigen.C"
#include "font_fix_hybrid_EOS.C"

