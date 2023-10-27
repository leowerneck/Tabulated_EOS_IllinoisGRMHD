//-----------------------------------------------------------------------------
// Compute the flux for advecting rho_star, tau (Font's energy variable),
//  and S_i .
//-----------------------------------------------------------------------------
static inline void mhdflux( const igm_eos_parameters eos,
                            const int i,
                            const int j,
                            const int k,
                            const int flux_dirn,
                            CCTK_REAL *restrict Ur,
                            CCTK_REAL *restrict Ul,
                            CCTK_REAL *restrict FACEVAL,
                            CCTK_REAL *restrict FACEVAL_LAPSE_PSI4,
                            CCTK_REAL &cmax,
                            CCTK_REAL &cmin,
                            CCTK_REAL &rho_star_flux,
                            CCTK_REAL &tau_flux,
                            CCTK_REAL &st_x_flux,
                            CCTK_REAL &st_y_flux,
                            CCTK_REAL &st_z_flux,
                            CCTK_REAL &Ye_star_flux,
                            CCTK_REAL &S_star_flux ) {

  CCTK_REAL psi4 = FACEVAL_LAPSE_PSI4[PSI4];
  CCTK_REAL psi6 = FACEVAL_LAPSE_PSI4[PSI4]*FACEVAL_LAPSE_PSI4[PSI2];
  CCTK_REAL psim4 = 1.0/(psi4);

  CCTK_REAL alpha_sqrt_gamma = FACEVAL_LAPSE_PSI4[LAPSE]*psi6;
  CCTK_REAL ONE_OVER_LAPSE = 1.0/FACEVAL_LAPSE_PSI4[LAPSE];
  CCTK_REAL ONE_OVER_LAPSE_SQUARED=SQR(ONE_OVER_LAPSE);

  // Compute the sound speed and enthalpy on the right face
  CCTK_REAL cs2_r=0,h_r=0;
  compute_cs2_and_enthalpy(eos,Ur,&cs2_r,&h_r);

  // Compute the sound speed and enthalpy on the left face
  CCTK_REAL cs2_l=0,h_l=0;
  compute_cs2_and_enthalpy(eos,Ul,&cs2_l,&h_l);

  //Compute face velocities
  // Begin by computing u0
  output_stats stats; stats.failure_checker=0;
  CCTK_REAL u0_r,u0_l;
  impose_speed_limit_output_u0(FACEVAL,Ur,psi4,ONE_OVER_LAPSE,stats,u0_r);
  impose_speed_limit_output_u0(FACEVAL,Ul,psi4,ONE_OVER_LAPSE,stats,u0_l);

  //Next compute b^{\mu}, the magnetic field measured in the comoving fluid frame:
  CCTK_REAL ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;
  /***********************************************************/
  /********** RIGHT FACE ************/
  // Note that smallbr[4] = b^a defined in Gammie's paper, on the right face.
  CCTK_REAL u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r;
  CCTK_REAL smallbr[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ur,u0_r,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4r,u_y_over_u0_psi4r,u_z_over_u0_psi4r,smallbr);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}] (UX=1,UY=2,UZ=3).
  CCTK_REAL U_LOWERr[4] = { 0.0, u_x_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4],
                            u_z_over_u0_psi4r*u0_r*FACEVAL_LAPSE_PSI4[PSI4] };
  /********** LEFT FACE ************/
  // Note that smallbl[4] = b^a defined in Gammie's paper, on the left face.
  CCTK_REAL u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l;
  CCTK_REAL smallbl[NUMVARS_SMALLB];
  // Compute b^{a}, b^2, and u_i over u^0
  compute_smallba_b2_and_u_i_over_u0_psi4(FACEVAL,FACEVAL_LAPSE_PSI4,Ul,u0_l,ONE_OVER_LAPSE_SQRT_4PI,
                                          u_x_over_u0_psi4l,u_y_over_u0_psi4l,u_z_over_u0_psi4l,smallbl);
  // Then compute u_xr,u_yr, and u_zr. We need to set the zeroth component so we can specify U_LOWER{r,l}[{UX,UY,UZ}]
  CCTK_REAL U_LOWERl[4] = { 0.0, u_x_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4], u_y_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4],
                            u_z_over_u0_psi4l*u0_l*FACEVAL_LAPSE_PSI4[PSI4] };
  /***********************************************************/

  // Compute v02 = v_A^2 + c_s^2*(1.0-v_A^2), where c_s = sound speed, and v_A = Alfven velocity
  CCTK_REAL v02r,v02l;
  // First right face
  compute_v02(eos,h_r,cs2_r,smallbr,Ur,v02r);
  // Then left face.
  compute_v02(eos,h_l,cs2_l,smallbl,Ul,v02l);

  int offset=flux_dirn-1;

  // Now that we have computed v02 = v_A^2 + c_s^2*(1.0-v_A^2), we can
  //   next compute c_+ and c_- using a simplified dispersion relation.
  //   Note that, in solving the simplified disp. relation, we overestimate
  //   c_+ and c_- by around a factor of 2, making the MHD evolution more
  //   diffusive (and potentially more *stable*) than it could be.
  CCTK_REAL cplusr,cminusr,cplusl,cminusl;
  find_cp_cm(cplusr,cminusr,v02r,u0_r,
             Ur[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);
  find_cp_cm(cplusl,cminusl,v02l,u0_l,
             Ul[VX+offset],ONE_OVER_LAPSE_SQUARED,FACEVAL[SHIFTX+offset],psim4,FACEVAL[GUPXX+offset]);

  // Then compute cmax, cmin. This is required for the HLL flux.
  CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
  CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));

  //*********************************************************************
  // density flux = \rho_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  CCTK_REAL rho_star_r = alpha_sqrt_gamma*Ur[RHOB]*u0_r;
  CCTK_REAL rho_star_l = alpha_sqrt_gamma*Ul[RHOB]*u0_l;
  CCTK_REAL Fr = rho_star_r*Ur[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ur[VX] -> Ur[VY]
  CCTK_REAL Fl = rho_star_l*Ul[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ul[VX] -> Ul[VY]

  // HLL step for rho_star:
  rho_star_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(rho_star_r-rho_star_l) )/(cmaxL + cminL);

  //*********************************************************************
  // Electron fraction flux = Ye_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  if( eos.is_Tabulated ) {
    CCTK_REAL Ye_star_r = rho_star_r * Ur[YEPRIM];
    CCTK_REAL Ye_star_l = rho_star_l * Ul[YEPRIM];
    Fr = Ye_star_r*Ur[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ur[VX] -> Ur[VY]
    Fl = Ye_star_l*Ul[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ul[VX] -> Ul[VY]

    // HLL step for Ye_star:
    Ye_star_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(Ye_star_r-Ye_star_l) )/(cmaxL + cminL);
  }

  //*********************************************************************
  // Entropy flux = S_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  if( eos.evolve_entropy ) {
    CCTK_REAL S_star_r = alpha_sqrt_gamma*Ur[ENTROPY]*u0_r;
    CCTK_REAL S_star_l = alpha_sqrt_gamma*Ul[ENTROPY]*u0_l;
    Fr = S_star_r*Ur[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ur[VX] -> Ur[VY]
    Fl = S_star_l*Ul[VX+offset]; // flux_dirn = 2, so offset = 1, implies Ul[VX] -> Ul[VY]

    // HLL step for S_star:
    S_star_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(S_star_r-S_star_l) )/(cmaxL + cminL);
  }

  //*********************************************************************
  // energy flux = \alpha^2 \sqrt{\gamma} T^{0m} - \rho_* v^m, where m is the current flux direction (the m index)
  //*********************************************************************
  // First compute some useful metric quantities:
  CCTK_REAL alpha_squared_sqrt_gamma = FACEVAL_LAPSE_PSI4[LAPSE]*alpha_sqrt_gamma;
  CCTK_REAL g4uptm =  ONE_OVER_LAPSE_SQUARED*FACEVAL[SHIFTX+offset];
  CCTK_REAL g4uptt = -ONE_OVER_LAPSE_SQUARED;
  /********** RIGHT FACE ************/
  // Compute a couple useful hydro quantities:
  CCTK_REAL rho0_h_plus_b2_r = (Ur[RHOB]*h_r + smallbr[SMALLB2]);
  CCTK_REAL P_plus_half_b2_r = (Ur[PRESSURE]+0.5*smallbr[SMALLB2]);
  // Then compute T^{0m} and the flux:
  CCTK_REAL TUP0m_r = rho0_h_plus_b2_r*SQR(u0_r)*Ur[VX+offset] + P_plus_half_b2_r*g4uptm - smallbr[SMALLBT]*smallbr[SMALLBX+offset];
  Fr = alpha_squared_sqrt_gamma * TUP0m_r - rho_star_r * Ur[VX+offset];
  // Finally compute tau
  CCTK_REAL TUP00_r = rho0_h_plus_b2_r*u0_r*u0_r + P_plus_half_b2_r*g4uptt - smallbr[SMALLBT]*smallbr[SMALLBT];
  CCTK_REAL tau_r = alpha_squared_sqrt_gamma * TUP00_r - rho_star_r;
  /********** LEFT FACE *************/
  // Compute a couple useful hydro quantities:
  CCTK_REAL rho0_h_plus_b2_l = (Ul[RHOB]*h_l + smallbl[SMALLB2]);
  CCTK_REAL P_plus_half_b2_l = (Ul[PRESSURE]+0.5*smallbl[SMALLB2]);
  // Then compute T^{0m} and the flux:
  CCTK_REAL TUP0m_l = rho0_h_plus_b2_l*SQR(u0_l)*Ul[VX+offset] + P_plus_half_b2_l*g4uptm - smallbl[SMALLBT]*smallbl[SMALLBX+offset];
  Fl = alpha_squared_sqrt_gamma * TUP0m_l - rho_star_l * Ul[VX+offset];
  // Finally compute tau
  CCTK_REAL TUP00_l = rho0_h_plus_b2_l*u0_l*u0_l + P_plus_half_b2_l*g4uptt - smallbl[SMALLBT]*smallbl[SMALLBT];
  CCTK_REAL tau_l = alpha_squared_sqrt_gamma * TUP00_l - rho_star_l;

  // HLL step for tau:
  tau_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(tau_r-tau_l) )/(cmaxL + cminL);


  //*********************************************************************
  // momentum flux = \alpha \sqrt{\gamma} T^m_j, where m is the current flux direction (the m index)
  //*********************************************************************
  // b_j = g_{ij} (b^i + b^t shift^i), g_{ij} = physical metric
  //CCTK_REAL sbtr=0,sbtl=0;
  CCTK_REAL smallb_lowerr[NUMVARS_SMALLB],smallb_lowerl[NUMVARS_SMALLB];
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbr,smallb_lowerr);
  lower_4vector_output_spatial_part(psi4,FACEVAL,smallbl,smallb_lowerl);

  /********** Flux for S_x **********/
  // [S_x flux] = \alpha \sqrt{\gamma} T^m_x, where m is the current flux direction (the m index)
  //    Again, offset = 0 for reconstruction in x direction, 1 for y, and 2 for z
  //    Note that kronecker_delta[flux_dirn][0] = { 1 if flux_dirn==1, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UX]
                          + P_plus_half_b2_r*kronecker_delta[flux_dirn][0] - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBX] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UX]
                          + P_plus_half_b2_l*kronecker_delta[flux_dirn][0] - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBX] );

  //        S_x =\alpha\sqrt{\gamma}( T^0_x )
  CCTK_REAL st_x_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UX] - smallbr[SMALLBT]*smallb_lowerr[SMALLBX] );
  CCTK_REAL st_x_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UX] - smallbl[SMALLBT]*smallb_lowerl[SMALLBX] );

  // HLL step for Sx:
  st_x_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_x_r-st_x_l) )/(cmaxL + cminL);

  /********** Flux for S_y **********/
  // [S_y flux] = \alpha \sqrt{\gamma} T^m_y, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][1] = { 1 if flux_dirn==2, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UY] + P_plus_half_b2_r*kronecker_delta[flux_dirn][1]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBY] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UY] + P_plus_half_b2_l*kronecker_delta[flux_dirn][1]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBY] );

  //        S_y =\alpha\sqrt{\gamma}( T^0_y )
  CCTK_REAL st_y_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UY] - smallbr[SMALLBT]*smallb_lowerr[SMALLBY] );
  CCTK_REAL st_y_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UY] - smallbl[SMALLBT]*smallb_lowerl[SMALLBY] );

  // HLL step for Sy:
  st_y_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_y_r-st_y_l) )/(cmaxL + cminL);

  /********** Flux for S_z **********/
  // [S_z flux] = \alpha \sqrt{\gamma} T^m_z, where m is the current flux direction (the m index)
  //    Again, offset = 1 for reconstruction in x direction, 2 for y, and 3 for z
  //    Note that kronecker_delta[flux_dirn][2] = { 1 if flux_dirn==3, 0 otherwise }.
  Fr = alpha_sqrt_gamma*( rho0_h_plus_b2_r*(u0_r*Ur[VX+offset])*U_LOWERr[UZ] + P_plus_half_b2_r*kronecker_delta[flux_dirn][2]
                          - smallbr[SMALLBX+offset]*smallb_lowerr[SMALLBZ] );
  Fl = alpha_sqrt_gamma*( rho0_h_plus_b2_l*(u0_l*Ul[VX+offset])*U_LOWERl[UZ] + P_plus_half_b2_l*kronecker_delta[flux_dirn][2]
                          - smallbl[SMALLBX+offset]*smallb_lowerl[SMALLBZ] );

  //        S_z =\alpha\sqrt{\gamma}( T^0_z )
  CCTK_REAL st_z_r = alpha_sqrt_gamma*( rho0_h_plus_b2_r*u0_r*U_LOWERr[UZ] - smallbr[SMALLBT]*smallb_lowerr[SMALLBZ] );
  CCTK_REAL st_z_l = alpha_sqrt_gamma*( rho0_h_plus_b2_l*u0_l*U_LOWERl[UZ] - smallbl[SMALLBT]*smallb_lowerl[SMALLBZ] );

  // HLL step for Sz:
  st_z_flux = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(st_z_r-st_z_l) )/(cmaxL + cminL);

  cmax = cmaxL;
  cmin = cminL;
}
