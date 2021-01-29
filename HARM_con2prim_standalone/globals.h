
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

 globals.h : 
 ---------------

    Purpose: To contain all the declarations for global variables,
              arrays, and functions.  This file is used to create
              "defs.h", a list of the definitions of these files.
              This file also contains all the definitions of the
              global structures.

              Users should not need to modify this file unless they
              are adding functionality.  Users may need to reference
              the definitions here in order to figure out how they are
              defined, e.g., what members a structure has.

              It is meant to be included in decs.h and therefore in
              most every source file.

******************************************************************************/

/*************************************************************************************/
/******   GLOBAL ARRAYS, VARIABLES, and FUNCTIONS ************************************/
/*************************************************************************************/

/**************************************************************************************
 These are global variables that store the (dynamic) values of these parameters that 
 were once MACROS.  Remove "_glob" from the name to reveal the original name of the 
 MACRO parameter 
***************************************************************************************/
extern int N1_glob;
extern int N2_glob;
extern int N3_glob;
extern int N1TOT_glob;
extern int N2TOT_glob;
extern int N3TOT_glob;
extern int NTOT_glob;
extern int N1E_glob;
extern int N2E_glob;
extern int N3E_glob;
extern int N_HIST_POINTS_glob;
extern int N_SURF_POINTS_glob;
extern int N_CONN_glob;
extern int CONN2_END1_glob;
extern int CONN2_END2_glob;
extern int CONN2_END3_glob;
extern int N1_CONN_glob;
extern int N2_CONN_glob;
extern int N3_CONN_glob;
extern int NM1_glob;
extern int NM2_glob;
extern int NM3_glob;
extern int N_GEOM_glob;
extern int MAX_NTOT_glob;
extern int N_COORD_glob;
extern int NM1E_glob;
extern int NM2E_glob;
extern int NM3E_glob;
extern int NM1_TOT_glob;
extern int NM2_TOT_glob;
extern int NM3_TOT_glob;
extern int N1_R_E_glob;
extern int N2_R_E_glob;
extern int N3_R_E_glob;
extern int N1_R_S_glob;
extern int N2_R_S_glob;
extern int N3_R_S_glob;
extern int DEGENERATE_glob;
extern int NM1_UP_BEG_glob;
extern int NM2_UP_BEG_glob;
extern int NM3_UP_BEG_glob;
extern int do_recon_dir[3];
extern int NCELLS_glob;


#if( MAKE_TIMERS ) 
extern double elapsed_times[N_TIMER_TYPES][TIMER_MICROSECONDS+1];
extern struct timeval beg_times[N_TIMER_TYPES];
extern struct timeval end_time;
extern int n_timer_interval;  /* number of steps between timer data dumps */
#endif


typedef unsigned long  int ulint; 
typedef unsigned short int usint; 

struct of_geom {
  double gcon[NDIM][NDIM] ;
  double gcov[NDIM][NDIM] ;
  double g ;
  double g_inv; 
  double alpha;
  double beta[NDIM-1];
  double ncon[NDIM];
} ;

struct of_state {
  double ucon[NDIM] ;
  double ucov[NDIM] ;
  double bcon[NDIM] ;
  double bcov[NDIM] ;
  double bsq;
  double p; 
  double eta;
  double ptot;
} ;

struct of_coord {
  int i,j,k,pos; 
  double x[NDIM];
  double xp[NDIM];
  double dx_dxp[NDIM][NDIM];
  double dxp_dx[NDIM][NDIM];
  double det_dx_dxp;
  double xcart[NDIM];
  double dxc_dxp[NDIM][NDIM];
  double r;
  double rhoflr, uuflr, tempflr, yeflr;
};

/* COORD_DIAGONAL3_DYN_RAD */
struct of_coord_diagonal3_dyn_rad { 
  double f0_r, fa_r, upsilon_r, Rin_of_t, dln_Rin_of_t_dt,t;
};

/* COORD_DIAGONAL4 */
struct of_coord_diagonal4 {
  int  n_exp_xp1;
  double h_xp2, delta, xp2_0, dth_min1, dth_min2;
  double r_jet, h_jet, exp_jet, r0_hyper, xp1_hyper;
  double c1, c2, t;
};

/* COORD_WARPED_SPHERICAL */
struct of_coord_warped_spherical { 
  double t;

  // warping parameters:
  double delta_x1, delta_x2, delta_x3, delta_x4 ;
  double delta_y1, delta_y2, delta_y3, delta_y4 ;
  double delta_z1, delta_z2 ;
  double a_x10, a_x20, a_z10 ;
  double h_x1, h_x2, h_x3, h_x4 ;
  double h_y1, h_y2, h_y3, h_y4 ;
  double h_z1, h_z2 ;
  double s_[3], br_[3] ;
  
  // Grid size
  double Rin_[3], Rout_[3] ;

  // time dependent and grid-point independent variables
  double ar_[3] ;
  double ybh_[3], xbh_[3], zbh_[3] ;
  double dxbh1_dt, dxbh2_dt, dybh1_dt, dybh2_dt, dybh3_dt;
  double d2xbh1_dt2, d2xbh2_dt2, d2ybh1_dt2, d2ybh2_dt2, d2ybh3_dt2 ;
  double rbh_[3] ;

  double square_per_x11 ;
  double square_per_x12 ;
  double square_per_x01 ;
  double square_per_x02 ;

  double square_int_per_zbh1;
  double square_int_per_z0  ;
  double square_int_per_z1  ;

  double square_int_per_xbh1;
  double square_int_per_xbh2;

  double square_int_per_x11 ;
  double square_int_per_x12 ;
  double square_int_per_x01 ;
  double square_int_per_x02 ;

  // time dependent and grid-point _dependent_ arrays
  double **square_xp1_3;
  double **square_xp1_4;

  double **dsquare_xp1_3;
  double **dsquare_xp1_4;

  double **square_per_xp2_1;
  double **square_per_xp2_2;
  double **square_per_xp3_1;
  double **square_per_xp3_2;
  double **square_per_xp3_3;
  double **square_per_xp3_4;

  double **square_int_per_xp2_1;
  double **square_int_per_xp2_2;
  double **square_int_per_xp3_1;
  double **square_int_per_xp3_2;
  double **square_int_per_xp3_3;
  double **square_int_per_xp3_4;

  double **dsquare_per_xp2_1;
  double **dsquare_per_xp2_2;

  double **dsquare_per_xp3_3;
  double **dsquare_per_xp3_4;

  double ***r_of_yt_;
  double ***drdy_of_yt_;
};

/* COORD_WARPED_CARTESIAN */
struct of_coord_warped_cartesian {
  double t;

  // warping parameters:
  double delta_x1, delta_x2, delta_x3, delta_x4 ;
  double delta_y1, delta_y2, delta_y3, delta_y4 ;
  double a_x10, a_x20, a_y10, a_y20 ;
  double h_x1, h_x2, h_x3, h_x4 ;
  double h_y1, h_y2, h_y3, h_y4 ;

  // Grid size
  double xmin, xmax, ymin, ymax ;
  
  // time dependent and grid-point independent variables
  double ybh_[3], xbh_[3], zbh_[3] ;
  double dxbh1_dt, dxbh2_dt, dybh1_dt, dybh2_dt ;
  //  double d2xbh1_dt2, d2xbh2_dt2, d2ybh1_dt2, d2ybh2_dt2 ;

  double square_per_x11 ;
  double square_per_x12 ;
  double square_per_x01 ;
  double square_per_x02 ;
  double square_int_per_xbh1;
  double square_int_per_xbh2;
  double square_int_per_x11 ;
  double square_int_per_x12 ;
  double square_int_per_x01 ;
  double square_int_per_x02 ;

  double square_per_y11 ;
  double square_per_y12 ;
  double square_per_y01 ;
  double square_per_y02 ;
  double square_int_per_ybh1;
  double square_int_per_ybh2;
  double square_int_per_y11 ;
  double square_int_per_y12 ;
  double square_int_per_y01 ;
  double square_int_per_y02 ;


  // time dependent and grid-point _dependent_ arrays

  double **square_per_xp2_1 ;
  double **square_per_xp2_2 ;
  double **square_per_xp2_3 ;
  double **square_per_xp2_4 ;
  double **dsquare_per_xp2_3 ;
  double **dsquare_per_xp2_4 ;

  double **square_int_per_xp2_1;
  double **square_int_per_xp2_2;

  double **square_per_xp1_1;
  double **square_per_xp1_2;
  double **square_per_xp1_3;
  double **square_per_xp1_4;

  double **dsquare_per_xp1_3 ;
  double **dsquare_per_xp1_4 ;

  double **square_int_per_xp1_1;
  double **square_int_per_xp1_2;
};

#if(COORD_TYPE_CHOICE==COORD_WARPED_SPHERICAL)
typedef struct of_coord_warped_spherical   of_coord_params;
#elif(COORD_TYPE_CHOICE==COORD_WARPED_CARTESIAN)
typedef struct of_coord_warped_cartesian   of_coord_params;
#elif(COORD_TYPE_CHOICE==COORD_DIAGONAL4)
typedef struct of_coord_diagonal4          of_coord_params;
#else
typedef struct of_coord_diagonal3_dyn_rad  of_coord_params;
#endif


/* Used for MPI masking function */
#if(USEMPI)
typedef struct { double val; int mask; } of_double_int; 
#endif

/* 123  MHD  functions per cell:     */
extern double    ****p  ;      /* space for primitive vars */
extern double    ****p_old;    /* previous time steps value of p[] */
extern double   *****p_L;      /* Left state for flux at i,j,k */
extern double   *****p_R;      /* slopes */
extern double   *****F;        /* fluxes */
extern double    ****ph;       /* half-step primitives */
extern double    ****U_gf[N0TOT];     /* Conserved variables */
extern double     ***p_gamma;  /* lorentz factor of latest primitive var's */
extern double    ****optical_depth;  /* optical depth for neutrinos */
extern double    ****opacity_nu;  /* opacity for neutrinos */
extern double   *****c_gf;     /* min/max wave speeds per face */
extern double   *****ucon_L;   /* ucon for entropy equation */
extern double   *****ucon_R;   /* ucon for entropy equation */
extern double     ***emf[NDIM-1];   /* Electromotive function for Constraint Transport method flux_ct */

extern double *p_vect, *dq_vect, *pL_vect, *pR_vect;


/* 148  GR functions per cell (connection not defined in ghosts) : */
extern struct of_geom *geom_arr[N0_GEOM]; 
extern double  ****conn;        /* connection coefficients */

#if(USE_PRESSURE_FLUX_FIX)
extern double   ****F2;        /* secondary flux */
#endif

/* 41*5 grid functions per cell: */
extern struct of_coord *coord_arr[N0_COORD];

extern int     ***pflag;          /* failure tracking array */
extern int     ***nfail[N_NFAIL]; /* numbers of different types of failures */

#if( USE_LOCAL_RECON_TYPE ) 
extern int    ****recon_type;     /* type of reconstruction for each cell */
#endif

extern usint ***evol_mask[N0_GEOM];  /* failure tracking array */

#if( RESCALE_REGULARIZE )
extern double *****regularize_gf,*****unregularize_gf;
#endif

#if( CALC_CURRENT ) 
extern double    ****jcon;     /* contravariant form of the current */
extern double  ******faraday;  /* gdet * (Faraday tensor) */
extern double    ****jcon2;    /* contravariant form of the current */
#endif

#if( USE_KINEMATIC_VISCOSITY )
struct of_viscosity {
  double eta_visc;
  double ucon[NDIM] ;
  double ucov[NDIM] ;
  double stress[NDIM][NDIM];  /* twice the Viscous stress tensor */
} ;
extern double               ****visc_source;    /* ultimate source term */
extern struct of_viscosity   ***visc_funcs[3];  /* covariant and contravariant 4-velocities */
#endif

#if( USE_COOLING_FUNCTION == 4 )
extern double  ***lum_thermal;  
extern double  ***dflux	     ;
extern double  ***dtau	     ;
extern int      **j_photo_bot;
extern int      **j_photo_top;
extern double   **tau_bot    ;
extern double   **tau_top    ;
extern double   **flux_top   ;
extern double   **flux_bot   ;
extern double   **flux_disk  ;
#if( USEMPI ) 
extern   double  *recv_buffer_top;
extern   double  *recv_buffer_bot;
extern   double  *send_buffer_top;
extern   double  *send_buffer_bot;
extern   double  *send_buffer_tot;
extern   double  *recv_buffer_tot;
extern   int    *irecv_buffer_top;
extern   int    *irecv_buffer_bot;
extern   int    *isend_buffer_top;
extern   int    *isend_buffer_bot;
extern   int    *isend_buffer_tot;
extern   int    *irecv_buffer_tot;
#endif
#endif


/* Statistical (STAT) gridfunctions : */
extern double *U_out[NP], *U_pre[NP], *dU_stat[NPH][N_STAT];
extern int stat_call_code;

/* STAT2 gridfunctions : */
extern double *FL_stat2[NDIM-1][NP],*FR_stat2[NDIM-1][NP];
extern double *UL_stat2[NDIM-1][NP],*UR_stat2[NDIM-1][NP];
extern double *ctop_stat2[NDIM-1], *S_stat2[2];

/* Grid parameters : */
extern double GridLength[NDIM]; /* Lengths of each side of simulation's spacetime volume    */
extern double startx[NDIM];     /* Coordinates of "earliest/leftmost" edge of domain volume */
extern double dx[NDIM], dV;     /* Step sizes and local cell's coordinate volume            */
extern double invdx[NDIM];      /* 1/dx                                                     */
extern double dt_old;           /* Previous step's time discretization                      */
extern double dtmin_glob;       /* global minimum light crossing time across a cell         */
extern double cour ;            /* Global Courant factor                                    */
extern double t,t_old;          /* Present, final, and previous times                       */
extern double dt_conn,half_inv_dt_conn; /* dt used in finite differencing metric for the connection */
extern int n_substep;           /* id number indicating which substep we are currently on   */
extern int n_U;                 /* time index of newest values of U_gf[]                    */

/* MPI-ish quantities: */
extern int    totalsize[NDIM];  /* Global dimensions; unique from N1-3 only if using MPI */
extern int    bc_pid[NDIM][2];  /* pid that local grid abuts against, or -NUM if physical BC  */
extern int    tag_send[2*3];    /* MPI tag id number for sending   data, indexed by face id   */
extern int    tag_recv[2*3];    /* MPI tag id number for receiving data, indexed by face id   */
extern int    myid;             /* MPI ID number for local grid     */
extern int    numprocs;         /* Total number of CPUs used in this run */
extern long   n_cells_glob;     /* Total number of physical cells in global domain (all of simulation) */
extern int    cpupos[NDIM];     /* Local grid's position (index) in each direction of global grid */
extern int    globalpos[NDIM];  /* Coordinates (in global grid indices) of local grid's 1st physical cell */
extern int    n_spatial_dims;   /* Number of independent spatial dimensions (e.g. 1D, 2D, 3D) */ 
extern int    ncpux[NDIM];      /* number of cpu's in each direction  */
extern int     printer_pid;     /* id of rank that reports messages; */
extern int          io_pid;     /* id of rank that gathers (if using GATHER_IO) and writes data to disk; */
extern int      master_pid;     /* id of rank that controls the entire job (usually equals "0"); */
extern int    special1_pid;     /* id of rank that performs special computations of type 1 */
extern int    ener_out_pid;     /* id of rank that write *.ener.out file   */
extern int    traj_out_pid;     /* id of rank that write *.ener.out file   */
extern int      timers_pid;     /* id of rank that is responsible for gathering and writing timing data */
extern int  out_pid[N_OUT_TYPES]; /* pid's to be responsible for that kind of IO type */

/* Variables controlling frequency of output types and IO arrays : */ 
extern int nstep;
extern int n_restart ;                /* Number of timesteps between restart dumps */
extern int rdump_cnt ;                /* Index number used to keep track of restart dumps */
extern int N_out[N_OUT_TYPES];        /* id numbers for next image/dump/rdump file */
extern double T_out[N_OUT_TYPES];     /* Time of next dump/hdf5/image */
extern double DT_out[N_OUT_TYPES] ;   /* Time intervals between various types of output data */ 
extern char DIR_out[N_OUT_TYPES][30]; /* Directories in which to store output files */
extern int N_hist_dump_tot;           /* Number of history dumps to do : */
extern int N_hist_dump;               /* Number of history dumps made so far */
extern double *hist_data[N_HISTORY][N_HIST_TYPES];
extern double *hist_data_out[N_HISTORY][N_HIST_TYPES];  /* array containing this timestep's history data */
extern double *hist_send;
extern double *hist_send2;
extern double **history_buffer;
extern double **history_buffer2;
extern double *surf_data[N_SURFACE][N_SURF_TYPES];
extern double *surf_send;
extern double *surf_recv;
  
extern double *f_sdf, *f_sdf2;
extern double *f_hdf[N_HDF_FUNCS];
extern double *U_avg[NP];

extern int n_r_bins, n_phi_bins;               /* numbers of radial and azimuthal bins for binning purposes  */
extern double *xp1_bins, *r_bins, *dr_bins, *inv_dr_bins, *phi_bins;   /* r,phi coordinates of lower edge of bins (i.e. vertex coordinates) */
extern double *r_bins_mid, *phi_bins_mid;   /* r,phi coordinates of middle of the bins (i.e. cell centered coordinates) */
extern double r0_bins, r_min_bins, r_max_bins; /* parameters defining the radial bin coordinates */
extern double xp1_min_bins, xp1_max_bins, dxp1_bins, inv_dxp1_bins;  /* parametesr for the uniform coordinates of the radial bins */
extern double phi_min_bins, phi_max_bins, dphi_bins, inv_dphi_bins;  /* parametesr for the uniform coordinates of the azimuthal bins */
extern double **bin_weights_r,**bin_weights_phi;
static int *ibins_beg, *ibins_end, *kbins_beg, *kbins_end;


/* Variables regarding failures */
extern short int failed;    /* Flag to indicate some type of failure that is usually algorithmic or unexpected */
extern short int sibling_failed;  /* Flag to indicate that a sibling process failed and not this process */
extern short int failure_exists ;  /* flag to indicate if an unphysical state has occurred */
extern int failure_type[N_PFLAG_TYPES];  /* to indicate which type of interpolation to do */
extern int treat_floor_as_failure ;  /* if not zero, will interpolate over regions with (rho,u < FLOOR) */ 
extern int using_restart; 
extern ulint nfailtot[N_NFAIL];
extern int boundary_mpi_pflag,boundary_phys_pflag; /* whether pflag was set within cells involved in mpi or physical boundary operations */

/* Variables that indicate position in spacetime: */
extern int  icurr, jcurr, kcurr, pcurr, ncurr;  /* Present cell & face indices */
extern int  n_beg, n_mid, n_end;  /* variables to store the beginning, middle (halfstep) and final time index for metric functions */
extern int  i_pscalar;   /* global index  for  passive scalar loop */
extern int  i_optdepth;   /* global index  for  optical depth  loop */

/* EOS variables */
extern double temp_guess;

/* Physics variables: */
extern double gam ;             /* Adiabatic index  */
extern double gam_m1_o_gam;
extern double *Katm;
extern double *nu_visc;     /* Kinematic viscosity */
extern double alpha_visc;         /* alpha ala "stress over pressure" parameter */
extern double E_tot;  /* E_tot =   \int_{domain}  T^t_t \sqrt{-g} dx^4  */
extern double M_tot;  /* M_tot =   \int_{domain}  rho   \sqrt{-g} dx^4  */
extern double L_tot;  /* L_tot =   \int_{domain}  T^t_t \sqrt{-g} dx^4  */
extern double *coolflux[NDIM];
extern double initial_bbh_separation;  /* initial separation of binary black hole (in units of m_bh1+m_bh2) */
extern double m_bh1;                   /* mass of BH1 (left BH)  */
extern double m_bh2;                   /* mass of BH2 (right BH) */
extern double m_bh_tot;                /* total mass of both BHs */
extern double M;                       /* mass scale, usually 1 or m_bh_tot */
extern double a_bh1;                   /* spin of BH1 */
extern double a_bh2;                   /* spin of BH2 */
extern double bh1_traj[N0_GEOM][NDIM];        /* trajectory of BH1 (left BH)  */
extern double bh2_traj[N0_GEOM][NDIM];        /* trajectory of BH2 (right BH)  */
extern double r_horizon1;             /* Radius of the horizon about BH1 */
extern double r_horizon2;             /* Radius of the horizon about BH2 */
extern double rsq_horizon1;           /* Square of Radius of the horizon about BH1 */
extern double rsq_horizon2;           /* Square of Radius of the horizon about BH2 */
extern int    MASS_TYPE;              /* Switch for mass value used by isotropic metric */
extern double t_shrink_bbh;           /* Time at which to start shrinking the binary */
extern double phi_0_bbh;              /* Phase of binary at t=0. */
#if( METRIC_DYNAMIC_TYPE_CHOICE!=METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )
extern double *nz_params; 
extern double nz_params_all[3][38];          
#endif // #if( METRIC_DYNAMIC_TYPE_CHOICE!=METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )

#if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND || METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND )
// IZ1 and IZ2 arrays for IZ functions of time and IZ constants:

extern double *iz1_params;
extern double *iz2_params;

extern double iz1_const;
extern double iz2_const;

// IZ1 and IZ2 tidal fields (functions of time only):
// Rank 2 Et:
extern double (*iz1_A2_Et)[4];
extern double (*iz2_A2_Et)[4];

// Rank 2 Etdot:
extern double (*iz1_A2_Etdot)[4];
extern double (*iz2_A2_Etdot)[4];

// Rank 2 Etdot0:
extern double (*iz1_A2_Etdot0)[4];
extern double (*iz2_A2_Etdot0)[4];

// Rank 3 Et:
extern double (*iz1_A3_Et)[4][4];
extern double (*iz2_A3_Et)[4][4];
// Rank 3 Ct:
extern double (*iz1_A3_Ct)[4][4];
extern double (*iz2_A3_Ct)[4][4];
// Rank 3 Ctdot:
extern double (*iz1_A3_Ctdot)[4][4];
extern double (*iz2_A3_Ctdot)[4][4];
// Rank 3 Ctdot0:
extern double (*iz1_A3_Ctdot0)[4][4];
extern double (*iz2_A3_Ctdot0)[4][4];
// Rank 4 Ct:
extern double (*iz1_A4_Ct)[4][4][4];
extern double (*iz2_A4_Ct)[4][4][4];

extern double iz1_params_all[3][200];
extern double iz2_params_all[3][200];
// Rank 2 Et:
extern double iz1_A2_Et_all[3][4][4];
extern double iz2_A2_Et_all[3][4][4];
// Rank 2 Etdot:
extern double iz1_A2_Etdot_all[3][4][4];
extern double iz2_A2_Etdot_all[3][4][4];
// Rank 2 Etdot0:
extern double iz1_A2_Etdot0_all[3][4][4];
extern double iz2_A2_Etdot0_all[3][4][4];
// Rank 3 Et:
extern double iz1_A3_Et_all[3][4][4][4];
extern double iz2_A3_Et_all[3][4][4][4];
// Rank 3 Ct:
extern double iz1_A3_Ct_all[3][4][4][4];
extern double iz2_A3_Ct_all[3][4][4][4];
// Rank 3 Ctdot:
extern double iz1_A3_Ctdot_all[3][4][4][4];
extern double iz2_A3_Ctdot_all[3][4][4][4];
// Rank 3 Ctdot0:
extern double iz1_A3_Ctdot0_all[3][4][4][4];
extern double iz2_A3_Ctdot0_all[3][4][4][4];
// Rank 4 Ct:
extern double iz1_A4_Ct_all[3][4][4][4][4];
extern double iz2_A4_Ct_all[3][4][4][4][4];


#endif //  METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND


#if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )

extern double chi_bh1;                   /* dimensionless BH1 spin (left BH)  */
extern double chi_bh2;                   /* dimensionless BH2 spin (right BH) */

// 1st order spin BBH (IZ+NZ) metric parameter structure:
  struct of_1st_order_spin_bbh_params {
    // Not explicitly in NZ:
    double t;       // Current simulation time
    // Not explicitly in NZ or FZ:
    double t_c;     // Coalescence time from initial shrinking time
    double phi;     // Current binary orbital phase (w.r.t. phi(t=0) = 0)
    double omega;   // Current binary orbital phase rate of change
    double r12;     // Current binary separation
    double r12dot;  // Current binary radial velocity (separation time derivative)
    // Buffer zone parameters (Not explicitly in NZ or FZ):
    double r1T0;    // Inner1-Near radius trans. func. parameter 
    double w1T0;    // Inner1-Near width trans. func. parameter 
    double r2T0;    // Inner2-Near radius trans. func. parameter 
    double w2T0;    // Inner2-Near width trans. func. parameter 
    double xNT0;    // Near-Inner radius trans. func. parameter 
    double wNT0;    // Near-Inner width trans. func. parameter 
    double lambda;  // Near-Far trans. func. parameter
    // Binary BH trajectory coordinates:
    double xi1x;
    double xi1y;
    double xi1z;
    double xi2x;
    double xi2y;
    double xi2z;
    // Binary BH velocity components:
    double v1x;
    double v1y;
    double v1z;
    double v2x;
    double v2y;
    double v2z;
    // Binary BH velocity magnitudes:
    double v1;
    double v2;
    // v1-v2 components:
    double v12x;
    double v12y;
    double v12z;
    // v1-v2 magnitude:
    double v12;
    // Dot product of v1 and v2:
    double v1v2;
    // Definitions used below: 
    // ni = (x - xi)/ |x - xi| 
    // is the unit direction to the field from source, and
    // ri = |x - xi| 
    // is the distance magnitude from the source to the field point.
    // n1-n2 components:
    double n12x;
    double n12y;
    double n12z;
    // Dot products:
    double n12v12; // (n1-n2).(v1-v2)
    double n12v1;  // (n1-n2).v1
    double n12v2;  // (n1-n2).v2
    // Spin definition: "Newtonian" spin S = c S_true = G m^2 chi
    // where S_true is the physical body spin which is 0.5PN order,
    // and chi is dimensionless spin parameter. The PN and PM metrics
    // on the NZ and FZ are usually written in terms of S.
    // Note that we use the dimensional Kerr spin parameter 
    // in the IZ metric: a = S/m.
    // Binary dimensionless spin components:
    double chi1x;
    double chi1y;
    double chi1z;
    double chi2x;
    double chi2y;
    double chi2z;
    // IZ1 and IZ2 tidal fields (functions of time only):
    double iz1_zRI[5];
    double iz2_zRI[5];
  };

typedef struct of_1st_order_spin_bbh_params of_bbh_params;
extern of_bbh_params *bbh_params;
extern of_bbh_params bbh_params_all[3];     // To account for 2nd order FDA operators 

  struct of_1st_order_spin_bbh_consts {
    double m1;       // BH 1 mass
    double m2;       // BH 2 mass
    double m;        // m1+m2
    double mu;       // m1*m2/m
    double eta;      // mu/m
    double delta_m;  // m1-m2
    double delta;    // delta_m/m
    double b;        // Initial binary separation
  // IZ1 and IZ2 constants (reserved constants):
    double iz1_const;
    double iz2_const;
  };

typedef struct of_1st_order_spin_bbh_consts of_bbh_consts;
extern of_bbh_consts bbh_consts;

#endif // #if( METRIC_DYNAMIC_TYPE_CHOICE==METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP )



/* grid dependent quantities : */
extern int   n_within_horizon; 
extern double a, asq, r_horizon, r_isco;
extern double x0_bh, y0_bh, z0_bh;    /* x,y,z coordinates of the black hole */
extern double R0, Rin, Rout;
extern double dx_global_min, dt_global_min;
extern double th_cutout, th_beg, th_end, th_length;
extern double h_slope, X1_0, X1_slope;                          /* COORD_DIAGONA and COORD_MIXED */
extern int    diag3_exponent;                                   /* COORD_DIAGONAL3 */
extern double diag3_factor;                                     /* COORD_DIAGONAL3 */
extern int    n_diag2_lines;                                    /* COORD_DIAGONAL2 */
extern double d0_diag2,*xi_diag2,*ci_diag2,*di_diag2,*si_diag2; /* COORD_DIAGONAL2 */

extern of_coord_params *coord_params;
extern of_coord_params  coord_params_all[3];
extern double x1_min, x1_max, x2_min, x2_max, x3_min, x3_max;
extern double *rscale[6];
extern double *inv_rscale[6];

/* phi-averaged metric parameters : */
extern int  n_phi_avg;   /* number of phi points used to perform the phi-average */
extern double d_phi_avg; /* resolution in phi */
extern double inv_n_phi_avg;  /*  1. / n_phi_avg */

#if(USE_NUMREL_DATA)
extern double t0_ML;

extern double *tj, *djdt;
extern int djdt_lines;
extern gsl_interp_accel *djdt_acc;
extern gsl_spline *djdt_spline;

extern double *tm, *dmdt;
extern int dmdt_lines;
extern gsl_interp_accel *dmdt_acc;
extern gsl_spline *dmdt_spline;

extern double *ttrack, *rtrack, *phitrack;
extern int track_lines;
extern gsl_interp_accel *rtrack_acc;
extern gsl_spline *rtrack_spline;
extern gsl_interp_accel *phitrack_acc;
extern gsl_spline *phitrack_spline;
#endif

#if(OUTPUT_MAX_VCHAR)
extern double ****max_vchar;
#endif


/*************************************************************************************/
/******   GLOBAL ROUTINES     ********************************************************/
/***** (note that all extern function declarations must be on one line ***************/
/*****  so that the makefile makes defs.h without mistakes) **************************/
/*************************************************************************************/

/* bounds.c : */ 
extern void bounds(double ****prim_arr, int skip_mpi);
extern void fix_flux( void ) ;

/* phys.c : */ 
extern void primtoflux(double *pr, struct of_state *q, int dir, struct of_geom *geom, double *flux) ;
extern void primtoU(double *pr, struct of_state *q, struct of_geom *geom, double *U);
extern void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) ;
extern void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd) ;
extern double bsq_calc(double *pr, struct of_geom *geom);
extern void get_state(double *pr, struct of_geom *geom, struct of_state *q);
extern void ucon_calc(double *pr, struct of_geom *geom, double *ucon);
extern int gamma_calc(double *pr, struct of_geom *geom, double *gamma);
extern double  my_sqrt_one_plus_x(double x); 
extern int gamma_calc2(double *pr, struct of_geom *geom, double *gamma);
extern void ucon2pr( double *pr, double ucon[NDIM], double gcon[][NDIM] );
extern void vchar(double *pr, struct of_state *q, struct of_geom *geom, int js, double *vmax, double *vmin);
extern void vchar_fast(struct of_state *q, struct of_geom *geom, int js, double *vmax, double *vmin);
extern void set_levi_civita(void);
extern void faraday_calc(int dim, struct of_state *q, double fcon[]);
extern void current_calc( void ) ;
extern double cooling_func_hr_disk( int i, int j, int k, double *ph );
extern void vel_circular_equatorial(double r, double vcon[NDIM]);
extern double omega_circular_equatorial(double r_loc);
extern double omega_circular_equatorial_general(double r_loc, double spin, double mass);
extern double cooling_func_isentropic_disk( int i,  int j, int k, double *ph, struct of_state *q);
extern double cooling_func_isentropic_disk_general( int i,  int j, int k, double *ph, struct of_state *q);
extern double target_temperature( double r_loc, double h_o_r );
extern void  calc_visc_stress_source( void );
extern void  set_visc_stress_functions(int nn, int ii, int jj, int kk, double *prim, struct of_state *q );
extern void  set_all_visc_stress_functions(int nn, double ****prim);
extern void  bound_visc_stress_functions(int nn, double ****prim);


/* metric.c : */
extern void lower(double *ucon, struct of_geom *geom, double *ucov);
// extern void raise(double *ucov, struct of_geom *geom, double *ucon);
extern void setutcon(double *vcon, double gcov[][NDIM]);
extern double risco_calc( int do_prograde );
extern double risco_calc_general( int do_prograde, double spin, double mass);
extern double rhorizon_calc(int pos_sign);
extern void gcov_func(struct of_coord *coords, double gcov[][NDIM]) ;
extern void bl_to_ks_cov(double *x, double blcov[], double kscov[] );
extern void bl_to_ks_con(double *x, double blcon[], double kscon[] );
extern void ks_to_bl_cov(double *x, double kscov[], double blcov[] );
extern void ks_to_bl_con(double *x, double kscon[], double blcon[] );
extern void bl_gcon_func( double *x, double gcon[][NDIM]);
extern void bl_gcov_func( double *x, double gcov[][NDIM]);
extern void bl_gcon_func_old( double *x, double gcon[][NDIM]);
extern void bl_gcov_func_old( double *x, double gcov[][NDIM]);
extern void bl_dxc_dxs_calc(double *x_cart, double *x_spher, double dxc_dxs[][NDIM] );
extern void bl_dxs_dxc_calc(double *x_spher, double *x_cart, double dxs_dxc[][NDIM] );
extern void ks_gcov_func(double *x, double gcov[][NDIM]);
extern void ks_gcon_func(double *x, double gcon[][NDIM]);
extern void ks_gdet_func(double *x, double *gdet );
extern void ks_cart_gcov_func(double *x, double gcov[][NDIM]);
extern void ks_cart_gcon_func(double *x, double gcon[][NDIM]);
extern void ks_cart_gdet_func(double *x, double *gdet );
extern void ks_cart_gcov_gcon_func(double *x, double gcov[][NDIM], double gcon[][NDIM]);
extern void ks_cart_to_ks_spher_pos(double *x_cart, double *x_spher);
extern void ks_dxc_dxs_calc(double *x_cart, double *x_spher, double dxc_dxs[][NDIM] );
extern void ks_dxs_dxc_calc(double *x_spher, double *x_cart, double dxs_dxc[][NDIM] );
extern void ks_cart_to_ks_spher_cov(double *x_cart, double *x_spher, double ks_cart_cov[], double ks_spher_cov[] );
extern void ks_cart_to_ks_spher_con(double *x_cart, double *x_spher, double ks_cart_con[], double ks_spher_con[] );
extern void ks_spher_to_ks_cart_pos(double *x_spher, double *x_cart);
extern void ks_spher_to_ks_cart_cov(double *x_spher, double *x_cart, double ks_spher_cov[], double ks_cart_cov[] );
extern void ks_spher_to_ks_cart_con(double *x_spher, double *x_cart, double ks_spher_con[], double ks_cart_con[] );
extern void ks_cart_conn_func(struct of_coord *coords, double ***connp );
extern void dxc_dxs_calc(double *x_cart, double *x_spher, double dxc_dxs[][NDIM] );
extern void dxs_dxc_calc(double *x_spher, double *x_cart, double dxs_dxc[][NDIM] );
extern void test_geom(void) ;
extern void get_special_geometry( struct of_coord *coords, struct of_geom *geom, int geom_type);
extern void set_general_conn( double t_input );
extern void set_excision_mask( usint  ***evol_mask_loc );
extern void init_general_metric(void);
extern void calc_all_geom(void);
extern void mink_cart_gcov_func(double *x, double gcov[][NDIM]);
extern void my_sincos( double th, double *sth, double *cth);

/* init.c : */
extern void init(void);
extern void init_base(void);


/* diag.c : */
extern void diag(int call_code);
extern void fail(int fail_type, int state);
extern double divb_cen_calc( int i, int j, int k ) ;
extern double divb_calc( int i, int j, int k ) ;
extern int current_time(void) ;
extern void trace_message(char *routine_name, int instance);


/* fixup.c : */ 
extern void fixup(double ****pv) ;
extern int check_floor(double *pv, struct of_coord *coords);
extern void fixup_floor(double *pv, struct of_coord *coords);
extern int check_gamma(double gamma);
extern void fixup_gamma(double *pv, double *gamma);
extern int check_Tmax(double *pv);
extern void fixup_Tmax(double *pv);
extern int check_entropy_eq(double rho, double uu, double bsq );
extern int fixup_entropy_eq(int i, int j, int k, double *prim_old, double *prim, struct of_geom *geom, struct of_geom *geom_old, double *gamma_out);
extern void set_Katm( void );


/* coord.c : */ 
extern void coord(int i, int j, int k, int loc, double *xp);
extern void all_x_of_xp( struct of_coord *coords) ;
extern void x_of_xp( double *x, double *xp ) ;
extern void dx_dxp_calc( double *x, double *xp, double dx_dxp[][NDIM] );
extern double det_dx_dxp_calc( double *x, double *xp );
extern double det_dx_dxp_calc2( double dx_dxp[][NDIM] );
extern void dx_dxp_dxp_calc( double *x, double *xp, double dx_dxp_dxp[][NDIM][NDIM] );
extern void transform_rank2cov(double *x, double *xp,  double gcov[][NDIM]);
extern void transform_rank2cov2( double dx_dxp[][NDIM], double vcov[][NDIM]);
extern void transform_rank2con( double *x, double *xp, double gcon[][NDIM]);
extern void transform_rank2con2( double dxp_dx[][NDIM], double vcon[][NDIM]);
extern void transform_rank1cov(double *x, double *xp,  double gcov[]);
extern void transform_rank1cov2(double dx_dxp[][NDIM], double vcov[]);
extern void transform_rank1con( double *x, double *xp, double gcon[]);
extern void transform_rank1con2( double dxp_dx[][NDIM], double vcon[]);
extern void transform_connection( double *x, double *xp, double conn[][NDIM][NDIM], double ***connp);
extern void transform_connection2( struct of_coord *coords, double conn[][NDIM][NDIM], double ***connp);
extern void inverse_transform_rank1con( double *x, double *xp, double vcon[]);
extern void inverse_transform_rank1con2(double dx_dxp[NDIM][NDIM],  double vcon[]);
extern void inverse_transform_rank1cov(double *x, double *xp,  double vcov[]);
extern void inverse_transform_rank1cov2(double dxp_dx[][NDIM], double vcov[]);
extern void inverse_transform_rank2con( double *x, double *xp, double vcon[][NDIM]);
extern void inverse_transform_rank2con2( double dx_dxp[][NDIM], double vcon[][NDIM]);
extern void inverse_transform_rank2cov(double *x, double *xp,  double vcov[][NDIM]);
extern void inverse_transform_rank2cov2(double dxp_dx[][NDIM], double vcov[][NDIM]);
extern void xcart_of_xspher(double *xcart, double *xspher, double dxc_dxs[][NDIM] );
extern void xcart_of_xspher_special(struct of_coord *coords );
extern void xcart_of_xspher_special2(struct of_coord *coords );
extern void xspher_of_xcart(double *xspher, double *xcart, double dxs_dxc[][NDIM] );
extern void xspher_of_xcart_only(double *xspher, double *xcart);
extern void xcart_of_xspher_only(double *xcart, double *xspher);
extern void xspher_of_xcart_special3(double *xspher, double *xcart, double *xp, double dxs_dxp[][NDIM] ) ;
extern void xcart_of_xspher_special3(double *xcart, double *xspher, double *xp, double dxc_dxp[][NDIM]) ;

extern double Rin_calc(void);
extern double find_min_dx(void);
extern double find_min_dt(void);
extern void check_cutout(void);
extern void check_coords( void );
extern void xcart_of_xspher(double *xcart, double *xspher, double dxc_dxs[][NDIM] );
extern void xcart_of_xspher_only(double *xcart, double *xspher);
extern void setup_diag2(int type);
extern void alloc_diag2(void);
extern void dealloc_diag2(void);
extern void regularize_rank1con( double *x, double vcon[NDIM-1]);
extern void unregularize_rank1con( double *x, double vcon[NDIM-1]);
extern void regularize_prim( double prim[NP], double reg[NDIM-1], double dxdxp[NDIM][NDIM] );
extern void unregularize_prim( double prim[NP], double unreg[NDIM-1], double dxpdx[NDIM][NDIM] );
extern void unregularize_prim_LR( double prim_L[NP], double prim_R[NP], double unreg[NDIM-1], double dxpdx[NDIM][NDIM] );
extern void calc_all_coord(int n, double t_now);
extern void  coord_of_r( double r, struct of_coord *coords ) ;
extern void  coord_of_xp( double *xp, struct of_coord *coords ) ;
extern void copy_coord(struct of_coord  *coords1, struct of_coord *coords2 ) ;
extern void update_coord_time_funcs(double t_now);
extern void set_bin_arrays();

extern void jac_xBL_to_xIZ( double Mass, double Spin, double xBL[NDIM], double dxIZ_dxBL[NDIM][NDIM] );
extern void jac_xIZ_to_xBL( double Mass, double Spin, double xIZ[NDIM], double dxBL_dxIZ[NDIM][NDIM]);
extern void jac_xNZ_to_xIZ( int BH,  double xNZ[NDIM], double dxIZ_dxNZ[NDIM][NDIM]);
extern void jac_xNZ_to_xIZ_r12const( int BH,  double xNZ[NDIM], double dxIZ_dxNZ[NDIM][NDIM]);
extern void coordsBL_of_coordsIZ( double Mass, double Spin, double xIZ[NDIM], double xBL[NDIM]);
extern void coordsIZ_of_coordsNZ( int BH, double xNZ[NDIM], double xIZ[NDIM]);
extern void transform_rank1con_BL2PN( int BH, struct of_coord *coords, double ucon[NDIM] );
extern void transform_rank1cov_BL2PN( int BH, struct of_coord *coords, double acov[NDIM] );

extern void get_stationary_coords(double *xcm, double *xref, double dxref_dxcm[NDIM][NDIM]);
/* restart.c : */
extern int restart_init(void);
extern void restart_write(void);


/* step_ch.c : */
extern void step_ch(void);
extern void step_ch_geom_only(void);
extern void step_ch_geom_only2(void);


/* harm_mpi.c */ 
extern void setup_mpi(int argc, char *argv[]);
extern void get_global_gfunc( int ip, double *gfunc ) ;
extern void set_mpi_grid(void);
extern void set_mpi_misc(void);
extern void myexit( int ret );
extern void get_global_ijk( int i, int j, int k, int  *i_glob, int *j_glob, int *k_glob) ;
extern int  get_global_index_pid( int pid, int i, int j, int k );
extern void write_chunk( int ivar,  double *f_chunk );
extern void mpi_global_min(double *minval);
extern void mpi_global_max(double *maxval);
extern void mpi_global_vmax(int n, double *maxval);
extern void mpi_global_sum(double *sumval);
extern void mpi_global_avg(double *avgval);
extern void mpi_global_imax(int *maxval);
extern void mpi_global_imax2(int *maxval1,int *maxval2);
extern void mpi_sum_to_head(double *v_send, double *v_recv, int n);
extern void mpi_ulong_sum_to_head(unsigned long int *v_send, unsigned long int *v_recv, int n);
extern void check_boundary_pflag( int i, int j, int k ) ;
extern void sync_val(double *val);
extern void sync_val_from_rank(double *val, int rank);
extern void sync_vect(double *val, int n);
extern void sync_vect_from_rank(double *val, int n, int rank );
extern void sync_int_vect_from_rank(int *val, int n, int rank );
extern void sync_int_from_rank(int *val, int rank);
extern void exit_status(void);


/* dump.c */
extern int c_to_fort_array( double *f_fort, double *f_c, int npts, int ndims, int *dims );

#if(USE_NUMREL_DATA)
extern void read_numrel_data( void );
extern void free_numrel_data( void );
extern void test_numrel_data( void );
#endif

/* decomp.c */
extern void set_global_domain_variables(void) ;

/* corona.c */
#if( USE_COOLING_FUNCTION == 4 )
extern void  setup_corona_cooling( double ****prim);
#endif
