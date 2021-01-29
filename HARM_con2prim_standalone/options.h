
/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

 options.h : 
 ---------------

    Purpose:  To contain all the compile-time  parameters that control the 
              properties of the simulation and how it is executed. 

              This is the only header file that a user is expected to
              modify when configuring a new run, unless they are
              adding functionality.

              It is meant to be included in decs.h and therefore most
              every source file.

******************************************************************************/

/* Note that the decision to use MPI is set by the variable "USEMPI" in the makefile.  */

/* Decide if we are using MPI-IO (!=0) or not (==0)  */
#define USE_MPI_IO      (0)

#define GATHER_IO       (1 && USEMPI)    /* whether or not to use 1-proc IO model, where data is gathered to master before writing */

/* Note that how it is setup now, history IO never uses MPI_IO even if this is unset: 
   (keep set so that one file is made) */
#define USE_MPI_IO_HIST (1)

/* Note that MPI-IO is not yet configured for OUT_SURFACE, so ignore this for now */
#define USE_MPI_IO_SURF (0)


/* Set this to be true if you want to produce tracing output : */
#define TRACE_CALLS (0) 


/* Set the number of seconds the run should last (a non-positive value allows it to run forever) */
//#define RUNTIME_SECONDS (172500)
#define RUNTIME_SECONDS (-1)

 /* Set to 1 if you want B^i = 0  and no evolution of B-field */
#define HYDRO_ONLY  (1)

/* Set to 1 if you want to avoid evolving the matter fields : */
#define NO_EVOLUTION (0)

/* Set to 1 if you want restrict the simulation to the spherical equator, need to setup the code in a particular way: 
   Typically used for debugging purposes, especially for debugging WARPED_SPHERICAL coordinates; */
#define EQUATORIAL_RUN (0)

/* Set to 1 if doing an equatorial run with a large cutout and an inspiralling binary. it will adjust the dt step size. */
#define LARGE_CUTOUT (0) 

#define CLEAN_MONOPOLES      (0)  /* Whether to clean div-B violations (aka magnetic monopoles) at the first time step:  */
#define INTERPOLATE_DATA     (0)  /* Whether to interpolate data (initial data or restart data) onto an arbitrary list of destination points (see interp_data.c) :  */
#define RUN_FROM_INTERP_DATA (0)  /* Whether to run from interpolated data */

/* Set to 1 if you want to output a file containing cartesian coordinates for IGM */
#define DUMP_COORDS_FOR_IGM (0)

/* Set to 1 if you want to run the simulation using the IGM Metric */
#define RUN_FROM_IGM_METRIC (0)

#define EVOLVE_ELECTRON_FRACTION  (1)   /* Whether or not to evolve the electron fraction and potentially use nuclear EOSs, opacities, emissivities... */

#define TABULATED_EOS_TABLES (1) /* Whether to use tabulates EOS tables or a gamma-law EOS*/

#define N_PASSIVE_SCALARS (2)    /* number of passive scalars you want to use */ 

#define N_OPTICAL_DEPTHS (4)    /* number of different optical depth functions to use, it should be twice the number of neutrino flavors included in the model */ 

#define RUN_TAG  "KDHARM0"   /* Set the name of the run: */

/*  Set this to true if the code should always abide to an array's index bounds, even when it is a multi-dimensional contiguous block of memory: */
#define USE_STRICT_ARRAY_BOUNDS  (1)

#define sincos my_sincos 

/*************************************************************************************/
/******   OPTIONS   and   CONSTANTS  (users should review this section)  *************/
/*************************************************************************************/

/* Choose the type of metric/coordinates we will be using:                               */
#define METRIC_MINK_CARTESIAN   (0)   /* Flatspace in cartesian x,y,z coordinates          */
#define METRIC_MINK_SPHERICAL   (1)   /* Flatspace in  spherical r,th,phi coordinates      */
#define METRIC_KS_SPHERICAL     (2)   /* Kerr-Schild in spherical r,th,phi coordinates     */
#define METRIC_KS_CARTESIAN     (3)   /* Kerr-Schild in cartesian x,y,z coordinates        */
#define METRIC_BL_SPHERICAL     (4)   /* Boyer-Lindquist in spherical r,th,phi coordinates */
#define METRIC_GENERAL_STATIC   (5)   /* Use a general metric defined externally but is NOT time-dependent  */
#define METRIC_GENERAL_DYNAMIC  (6)   /* Use a general time dependent metric               */
#define METRIC_GENERAL_PHI_AVG  (7)   /* Use a general time dependent metric but phi-avgerage it as if its a time-average */
#define METRIC_GENERAL_PHI_AVG2 (8)   /* Like METRIC_GENERAL_PHI_AVG2, but it is not transformed to xp coordinates; used for initial data only  */

#define METRIC_TYPE_CHOICE (METRIC_KS_SPHERICAL)


/* Choose the type of dynamical metric:              */
#define METRIC_DYNAMIC_KS_SPHERICAL                (0)  /* Kerr-Schild in spherical r,th,phi coordinates (for testing purposes) */
#define METRIC_DYNAMIC_MINK_CARTESIAN              (1)  /* For testing purposes with DYNAMIC_COORDINATES */
#define METRIC_DYNAMIC_MINK_SPHERICAL              (2)  /* For testing purposes with DYNAMIC_COORDINATES */
#define METRIC_DYNAMIC_35PN_NZ                     (3)  /* 3.5PN order near zone approximate to binary black hole spacetime */
#define METRIC_DYNAMIC_FULLPN_NZ                   (4)  /* "Full" PN order near zone approximate to binary black hole spacetime */
#define METRIC_DYNAMIC_FULLPN_NZ_FAST              (5)  /* "Full" PN order near zone approximate to binary black hole spacetime -- optimized further by Bruno Mundim */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ                (6)  /* "Full" PN order near zone and Inner Zone matched together */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ             (7)  /* "Full" PN order near zone, Inner Zone and Far Zone matched together */
#define METRIC_DYNAMIC_FULLPN_NR_FZ                (8)  /* "Full" PN order Far Zone along with Numerical Relativity merger and ringdown data */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_NEW         (9)  /* new "full" PN order near zone */
#define METRIC_DYNAMIC_FULLPN_NZ_DROP             (10)  /* "Full" PN order near zone with switches to drop the order */
#define METRIC_DYNAMIC_SIMPLE_NEWTONIAN           (11)  /* simple newtonian metric */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ_2ND           (12)  /* "Full" PN order near zone and Inner Zone matched together to 2nd order. IZ highly optimized now. */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_2ND        (13)  /* "Full" PN order near zone, Inner Zone and Far Zone matched together to 2nd order. IZ highly optimized now. */
#define METRIC_DYNAMIC_FULLPN_NZ_IZ_FZ_ASPIN_DROP (14)  /* "Full" PN order near zone, Inner Zone and Far Zone matched together to 1ST order. Aligned spins only. */

#define METRIC_DYNAMIC_TYPE_CHOICE (METRIC_DYNAMIC_KS_SPHERICAL)

/* Evolution and metric is in a comoving frame */
#define COMOVING_METRIC (0)    /* Transforms the metric calculated using the above options to a moving frame --Must be used with METRIC_GENERAL_DYNAMIC */

/* Choose boundary conditions ala global types (eg. cartesian, spherical, cylindrical)   */
#define BC_CARTESIAN_OUTFLOW     (0)      /* Outflow  on all boundaries             */
#define BC_CARTESIAN_PERIODIC    (1)      /* Periodic on all boundaries             */ 
#define BC_CARTESIAN_REFLECT     (2)      /* Periodic on all boundaries             */ 
#define BC_SPHERICAL_OUTFLOW     (3)      /* Outflow in X1, constant in X2, Per. in X3 */
#define BC_CYLINDRICAL_OUTFLOW   (4)      /* Outflow in X1,X2 , Periodic in X3      */
#define BC_SPHERICAL_CUTOUT      (5)      /* Outflow in X1, refl. in cutout along X2, Per. in X3 */
#define BC_SPHERICAL_CONSINTERP  (6)      /* Like OUTFLOW except X2 boundary copy uses U[]  */ 
#define BC_SPHERICAL_OUTFLOW2    (7)      /* Outflow in X1, axi-sym in X2, Per. in X3 */
#define BC_SPHERICAL_CUTOUT2     (8)      /* Outflow in X1, axi-sym in X2, Per. in X3 */
#define BC_SPHERICAL_CUTOUT3     (9)      /* Outflow in X1, axi-sym in X2, Per. in X3 */
#define BC_SPHERICAL_OUTFLOW3    (10)     /* Outflow in X1, reflective axi-sym in X2, Per. in X3 */
#define BC_SPHERICAL_OUTFLOW4    (11)     /* Outflow in X1, reflective axi-sym in X2, Per. in X3, no imposement of outflow at radial boundaries */
#define BC_SPHERICAL_OUTFLOW5    (12)     /* Outflow in X1, axi-sym in X2, Per. in X3, no imposement of outflow at radial boundaries */
#define BC_SPHERICAL_OUTFLOW6    (13)     /* Outflow in X1, reflective axi-sym in X2, Per. in X3, no imposement of outflow at inner radial boundaries, use inflow_check() at outer radial boundary*/
#define BC_SPHERICAL_OUTFLOW7    (14)     /* Outflow in X1, reflective axi-sym in X2, Per. in X3, no imposement of outflow at inner radial boundaries, use inflow_check2() at outer radial boundary*/
#define BC_STATIC_ALL            (15)     /* Boundary data kept constant */
#define BC_SPHERICAL_STATIC      (16)     /* Only radial BC's kept constant, reflective anti-sym in X2, Per. in X3, like BC_SPHERICAL_OUTFLOW3 */
#define BC_USER_SUPPLIED         (17)     /* Put your own routine in bounds.c       */ 
#define BC_SPHERICAL_FULLSPHERE  (18)     /* Spherical grid covers the spherical domain, sharing ghostzone data across axis and origin, uses special origin cell, outer radial outflow BC */
#define BC_SPHERICAL_FULLSPHERE2 (19)     /* Spherical grid covers the spherical domain, sharing ghostzone data across axis and origin, no origin cell, send data across origin, outer radial outflow BC */
#define BC_SPHERICAL_FULLSPHERE3 (20)     /* Spherical grid covers the spherical domain, sharing ghostzone data across axis and origin though setting v=0 at origin, no origin cell, only one outflow condition on the outer radial boundary */
#define BC_SPHERICAL_FULLSPHERE4 (21)     /* Spherical grid covers the spherical domain, sharing ghostzone data across axis and phi boundary, no origin cell, outflow conditions on both radial boundaries */
#define BC_SPHERICAL_OUTFLOW8    (22)     /* Outflow in X1, "EQUATORIAL_RUN=1" boundaries, i.g. copies in X2, Per. in X3, no imposement of outflow at inner radial boundary */
#define BC_SPHERICAL_OUTFLOW9    (23)     /* Outflow in X1, "EQUATORIAL_RUN=1" boundaries, i.g. copies in X2, Per. in X3, IMPOSEMENT of outflow at inner radial boundary */
#define BC_TYPE_CHOICE           (BC_SPHERICAL_FULLSPHERE4)


/* Choose how the "standard/BL/KS/etc" coordinates are defined w.r.t. the numerical coordinates */
#define COORD_IDENTITY           (0)      /* No difference, numerical = standard      */
#define COORD_DIAGONAL           (1)      /* "Classic" HARM or "old" MKS  coordinates */
#define COORD_MIXED              (2)      /* "New MKS" coordinates to resolve jets    */ 
#define COORD_DIAGONAL2          (3)      /* Like DIAGONAL except allows any sequence of th discretizations */
#define COORD_DIAGONAL3          (4)      /* Like DIAGONAL except th is a polynomial function of x2 */
#define COORD_DIAGONAL3_DYN_RAD  (5)      /* Dynamic version of COORD_DIAGONAL3 -- radial boundary moves */
#define COORD_WARPED_SPHERICAL   (6)      /* Dynamic Warped spherical coordinate system  (for binary black holes) */
#define COORD_WARPED_CARTESIAN   (7)      /* Dynamic Warped cartesian coordinate system */
#define COORD_DIAGONAL3_CUTOUT   (8)      /* The version of DIAGONAL3 to use when needing a sizeable polar axis cutout */
#define COORD_DIAGONAL4          (9)      /* Similar to DIAGONAL3 but WARPED_SPHERICAL's theta grid, and hyper-exponential radial grid */

#define COORD_TYPE_CHOICE  (COORD_DIAGONAL4)


/* whether to make sure that the coordinates make sense near for the ghost cells near the axis */
#define ENFORCE_GHOST_COORDINATE_CONDITIONS_AT_AXIS (0)
  

/* Reconstruction Method for "left/right" states of each cell's Riemann problem  */
/*  (note that non-linear methods should be >= RECON_PPM_MM ) */
#define RECON_MINMOD  (0)    /* 1st-order recon: mindmod slope limiter           */
#define RECON_VANL    (1)    /* 1st-order recon: van-leer slope limiter          */
#define RECON_MC      (2)    /* 1st-order recon: monotized central slope limiter */
#define RECON_PPM_MM  (3)    /* 2nd-order recon: piecewise parabolic method      */
#define RECON_PPM_VL  (4)    /* 2nd-order recon: piecewise parabolic method      */
#define RECON_PPM_MC  (5)    /* 2nd-order recon: piecewise parabolic method      */

#define RECON_TYPE_CHOICE   (RECON_PPM_MC)


/* Choose 1 so that we can change the reconstruction method cell-to-cell  */
#define USE_LOCAL_RECON_TYPE (0) 


/* Decide to use parabolic reconstruction of EMFs for Flux CT algorithm : */
#define USE_FLUX_CT_PARA (1)


/* Numerical Flux method */
#define FLUX_LF      (0)
#define FLUX_HLLE    (1)

#define FLUX_TYPE_CHOICE  (FLUX_LF)


/* Time step determining method: */
#define TIMESTEP_ORIGINAL            (0)    /*  original method */
#define TIMESTEP_MIN_CROSSING_TIME   (1)    /*  use minimum wave crossing time anywhere on grid */

#define TIMESTEP_METHOD  (TIMESTEP_ORIGINAL)


 /* Whether or not to use light speed or MHD characteristic speed to calculate the global minimum light crossing time across a cell : */ 
#define USE_LIGHT_SPEED  (0) 


/* whether to use the cooling function in the source terms:
   0  ->  no cooling fuction
   1  ->  cool toward target temperature to maintain constant h/r 
   2  ->  cool toward constant target entropy using local Keplerian period in cooling rate
   3  ->  cool toward constant target entropy using u^\phi/u^t in cooling rate
   4  ->  like 1 but also cool corona at an approximate inv. Compton rate 
 */
#define USE_COOLING_FUNCTION (0)   

 /* Whether or not to separate pressure term from fluxes to cancel source term that becomes a problem near the origin */ 
#define USE_PRESSURE_FLUX_FIX  (0) 


 /* Whether to use kinematic viscosity :  = 1 (general)  ,   = 2 (only T^r_\phi viscosity)  */ 
#define USE_KINEMATIC_VISCOSITY (0) 


#define FIXUP_TREE  (3)   /* Selects the sequence of fixups, see ./info/fixups */
                          /* 3 -> like 2 but never changes the B-field components */


#define ALLOW_NEGATIVE_DENSITIES  (0)   /* Allow rho,u to be negative  */

/* Option to use the entropy evolution equation instead of total energy equation */ 
#define USE_ENTROPY_EQ   (1)       /* Whether to ever use entropy equation */
#define BETA_MIN         (1.e-2)   /* use EE fix when u < BETA_MIN * bsq */


#define MASS_BH          (3.0)  // mass of the BH, for units

#define KEEP_CONSERVED_VARS (0)  /* Whether to keep U[] around for next timestep */  //--BROKEN CODE WILL SEGFAULT

/* Option to use constrain uu so that uu/rho is never larger than TEMPERATURE_MAX : */
#define USE_TEMPERATURE_MAX  (0 )  
#define TEMPERATURE_MAX      (100.)  

/* Use a simplified equation of state for primitive variable calculation 
     Note: the choices in EOS are specified in u2p_defs.h and include only "isothermal" and "isentropic" */
#define USE_SIMPLE_EOS (0)

/* Floor types: (see fixup1zone() for details) */
#define FLOOR_STATIC     (0)    /* Floor profile that has most often been used with HARM2d */
#define FLOOR_FLAT       (1)    /* Floor profile that has been used by Hawley et al.       */
#define FLOOR_TEMP       (2)    /* Floor using the minimum temperature (currently set to TEMPMIN=1e7K), RHOMIN and ye=0.1*/

#if (TABULATED_EOS_TABLES)
  #define FLOOR_TYPE_CHOICE (FLOOR_TEMP)
#else
  #define FLOOR_TYPE_CHOICE (FLOOR_STATIC)  
#endif

/* Depending on different floor types, set the magnitude and scaling of the floor : */
#if( FLOOR_TYPE_CHOICE == FLOOR_STATIC ) 
#define RHOMIN      (2.e-10)    /* Floor on rest-mass density */
#define UUMIN       (2.e-12)    /* Floor on internal energy density */
//#define RHOMIN      (5.e-11)   /* Floor on rest-mass density */
//#define UUMIN       (6.e-16)   /* Floor on internal energy density */
#define RHOMINLIMIT (1.e-30)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define UUMINLIMIT  (1.e-30)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define RHOPOWER    (-1.5)     /* Power in power-law that rho in floor follows  */
#define UUPOWER     (-2.5)     /* Power in power-law that uu in floor follows  */
#define FLOOR_DIM   (COORD_RADIUS_DIM)  /* Dimensionality of floor grid array */

#elif( FLOOR_TYPE_CHOICE == FLOOR_FLAT ) 
//#define RHOMIN      (5.e-11)   /* Floor on rest-mass density */
//#define UUMIN       (6.e-18)   /* Floor on internal energy density */
#define RHOMIN      (1.e-8)   /* Floor on rest-mass density */
#define UUMIN       (1.e-10)   /* Floor on internal energy density */
#define RHOMINLIMIT (1.e-22)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define UUMINLIMIT  (1.e-22)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define RHOPOWER    (0.)       /* Power in power-law that rho in floor follows  */
#define UUPOWER     (0.)       /* Power in power-law that uu in floor follows  */
#define FLOOR_DIM   (0)        /* Dimensionality of floor grid array */
#elif( FLOOR_TYPE_CHOICE == FLOOR_TEMP ) 
//#define RHOMIN      (5.e-11)   /* Floor on rest-mass density */
//#define UUMIN       (6.e-18)   /* Floor on internal energy density */
#define RHOMIN      (5.e-14)   /* Floor on rest-mass density */
#define UUMIN       (-2.5e-03)   /* Floor on internal energy density it is the value using the current table of TEMPMIN=1e-3MeV */
#define TEMPMIN     (1.1e-1)     /* Floor on Temperature, in MeV*/
#define YEMIN       (0.005)
#define RHOMINLIMIT (1.e-22)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define UUMINLIMIT  (1.e-22)   /* Used for r-dep. floor, abs. min. of RHOMIN */
#define RHOPOWER    (-1.5)       /* Power in power-law that rho in floor follows  */
#define UUPOWER     (-2.5)       /* Power in power-law that uu in floor follows  */
#define TEMPPOWER   (-1.)       /* Power in power-law that uu in floor follows  */
#define FLOOR_DIM   (0)        /* Dimensionality of floor grid array */
#endif

#define FIX_CUTOUT_PRESSURE (0)  /* Whether to smooth/interpolate bad pressures along the cutout boundary */
                                 /* Only can be used with BC_SPHERICAL_CUTOUT (see fixup.c) */

#define FIX_CUTOUT_GAMMA    (0)  /* Whether to smooth/interpolate bad pressures along the cutout boundary */
                                 /* Only can be used with BC_SPHERICAL_CUTOUT (see fixup.c) */

#define USE_GAMMA_CEILING_X1DN_BC  (0)  /* Whether to limit gamma in ghost zones at X1DN boundary */


#define GAMMAMAX            (50.)   /* Ceiling for the lorentz factor */
#define GAMMAMAX2           (1.e4)  /* To indicate inversion failure */
#define RESCALE_B           (0)     /* Whether to reconstruct using (det(-g) * B^i) */
#define RESCALE_R           (1)     /* Whether to rescale primitives in r before recon */
#define RESCALE_REGULARIZE  (0)     /* Whether to regularize or remove geometric factors from prim. vectors (set to 2 to use in bounds())*/

#define COORDSINGFIX (1)  /* Set to 1 to use coord. fixes assoc. w/ roundoff errors at axis */
#define SINGSMALL (1.E-16)  /* closest thing to zero that theta gets */

#define FAST_AND_FURIOUS    (1)    /* Whether to skip some MPI communications used to check consistency and MPI alive states that may or may not be needed */

/* Excision Mask Macros: */
#define USE_MASK     (0)      /* Whether to use an mask that identifies an excision region */
#define BUFFER_WIDTH (3)      /* what number of cells is the buffer wide?   */

#define EXCISE_HORIZON (0)    /* The default is to excise the horizon when using mask */
#define EXCISE_IZ      (1)    /* Excise the Inner-zone instead */
#define EXCISE_ISCO    (2)    /* Excise a region bounded by ISCO of an isolated BH */
#define EXCISION_TYPE_CHOICE (EXCISE_HORIZON)

#if( USE_MASK )
#define EXTEND_CENTRAL_CUTOUT (0) /* Only implemented for BC_SPHERICAL_OUTFLOW6 currently */
#endif


 /**************** OUTPUT OPTIONS: ***********************/ 
#define MAKE_HISTORY   (0)  /* Whether to generate history    files */
#define MAKE_ASCII     (0)  /* Whether to generate ascii dump files */
#define MAKE_SDF       (0)  /* Whether to generate sdf        files */
#define MAKE_HDF5      (1)  /* Whether to generate hdf5       files (affects dump_history file too) */
#define MAKE_IMAGE     (0)  /* Whether to generate rasterized images of data slices */
#define MAKE_GDUMP     (0)  /* Whether to generate coordinate/metric functions */
#define MAKE_STAT      (0)  /* Whether to generate statistical information files (hdf only) */
#define MAKE_STAT2     (0)  /* Whether to dump left/right fluxes, U, P, sources, ctop (hdf only) */
#define MAKE_RADFLUX   (0)  /* Whether to make separate hdf5 file for radiative fluxes */
#define MAKE_MIN_DT    (0)  /* Whether calculate and print the shortest light crossing time on the grid, writes to stdout */
#define MAKE_SURFACE   (0)  /* Whether calculate and print the shortest light crossing time on the grid, writes to stdout */
#define MAKE_TIMERS    (0)  /* Whether to calculate and write   */
#define MAKE_EHT_DATA  (0)  /* Whether to calculate data suggested for the EHT code comparison project    */


#define OUTPUT_CONSERVED_VARS (0) /* Whether to store and dump the conserved variables */ 
#define OUTPUT_COORDS         (1) /* Whether to dump coordinates in 3d dumps */
#define OUTPUT_CONN           (0) /* Whether to dump connection coefficients in 3d dumps */
#define OUTPUT_MAX_VCHAR      (1) /* Whether to local fastest characteristic speeds used for the fluxes and timestep determination */
#define DUMP_ALL_STAT (0)  /* Whether to calculate and dump dU_stat[], U_pre, U_fin, U_avg */


#define CHECK_INVERSION (0)       /* Whether or not to measure the inversion error in U per zone-cycle */
#define INVERSION_THRESH (1.e-4)  /* Threshold above which to report rel. err. in U  */

#define CHECK_PRIM_CHANGE (0)     /* Whether to test for large changes in prim. var's thru inversion */
#define PRIM_CHANGE_THRESH (0.99)   

#define CALC_CURRENT (0)  /* Choose to calculate the current or not (uses pface[]) */


/* Method of calculating connection */ 
#define CONN_METHOD_ORIG       (0)
#define CONN_METHOD_2ND_ORDER  (1)
#define CONN_METHOD_4TH_ORDER  (2)

#define CONN_METHOD_CHOICE (CONN_METHOD_4TH_ORDER)


/* For dynamic spacetimes, set whether (1) or not (0) to start shrinking the binary at a future time and not immediately: */
#define SET_TSHRINK (0)   // use t_shrink_bbh parameter

/* Read in time dependent single black hole data from numerical relativity calculations for post-merger simulations:  */ 
#define USE_NUMREL_DATA   (0) 


