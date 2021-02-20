// Thorn      : IllinoisGRMHD
// File       : EOS_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains function prototypes for the
//              EOS functions available in IllinoisGRMHD

#ifndef IGM_EOS_HEADERS_HH
#define IGM_EOS_HEADERS_HH

// These are useful for us, though not defined
// explicitly by EOS_Omni
static const CCTK_INT table_key_pressure = 0;
static const CCTK_INT table_key_epsilon  = 1;
static const CCTK_INT table_key_entropy  = 2;

// These are also very useful
static const CCTK_INT have_eps  = 0;
static const CCTK_INT have_temp = 1;
static const CCTK_INT have_ent  = 2;
static const CCTK_INT have_prs  = 3;

//------------------- EOS struct -----------------

// This struct defines all the necessary quantities
// that IllinoisGRMHD needs in order to use the
// EOS functions correctly

// Maxmimum number of polytropic EOS pieces
#define MAX_EOS_PARAMS (10)

typedef struct _igm_eos_parameters_ {

  //-------------- General parameters --------------
  // EOS key (selects Hybrid or Tabulated EOS)
  CCTK_INT key;
  // This is useful so we know which EOS type we are using
  bool is_Hybrid,is_Tabulated;
  // con2prim routine key (selects the con2prim routine)
  CCTK_INT c2p_routine;
  // con2prim backup routines keys
  CCTK_INT c2p_backup[3];
  // Which variable to reconstruct in PPM
  CCTK_INT PPM_reconstructed_var;
  // Whether or not to evolve the entropy
  bool evolve_T,evolve_entropy;
  // Baryonic density parameters
  CCTK_REAL rho_atm, rho_min, rho_max;
  // Atmospheric tau
  CCTK_REAL tau_atm;
  // Maximum Lorentz factor
  CCTK_REAL W_max,inv_W_max_squared;
  //------------------------------------------------

  //----------- Hybrid Equation of State -----------
  // Number of polytrope pieces
  CCTK_INT neos;
  // Density values which devide the polytropic pieces
  CCTK_REAL rho_ppoly_tab[MAX_EOS_PARAMS-1];
  // Polytropic indices
  CCTK_REAL Gamma_ppoly_tab[MAX_EOS_PARAMS];
  // Adiabatic constants
  CCTK_REAL K_ppoly_tab[MAX_EOS_PARAMS];
  // Integration constants for specific internal energy
  CCTK_REAL eps_integ_const[MAX_EOS_PARAMS];
  // Thermal adiabatic index
  CCTK_REAL Gamma_th;
  //------------------------------------------------

  //---------- Tabulated Equation of State ---------
  // Electron fraction parameters
  CCTK_REAL Ye_atm , Ye_min , Ye_max ;
  // Temperature parameters
  CCTK_REAL T_atm  , T_min  , T_max  ;
  // Pressure parameters
  CCTK_REAL P_atm  , P_min  , P_max  ;
  // Specific internal energy parameters
  CCTK_REAL eps_atm, eps_min, eps_max;
  // Entropy parameters
  CCTK_REAL S_atm  , S_min  , S_max  ;
  // Root-finding precision for table inversions
  CCTK_REAL root_finding_precision;
  // This threshold is used by the Palenzuela con2prim routine
  CCTK_REAL depsdT_threshold;
  //------------------------------------------------

} igm_eos_parameters;
//------------------------------------------------

void initialize_igm_eos_parameters_from_input( const CCTK_INT* igm_eos_key, const CCTK_REAL cctk_time, igm_eos_parameters &eos );
void apply_floors_and_ceilings_to_prims__recompute_prims( const igm_eos_parameters eos, const CCTK_REAL *restrict METRIC_LAP_PSI4, CCTK_REAL *restrict PRIMS );

//----------- Hybrid Equation of State -----------
void print_EOS_Hybrid( igm_eos_parameters eos );
void setup_K_ppoly_tab__and__eps_integ_consts( igm_eos_parameters &eos );
void initialize_Hybrid_EOS_parameters_from_input( igm_eos_parameters &eos );
int find_polytropic_K_and_Gamma_index( igm_eos_parameters eos, CCTK_REAL rho_in );
int find_polytropic_K_and_Gamma_index_from_P( const igm_eos_parameters eos, const CCTK_REAL P_in );
void compute_P_cold__eps_cold(igm_eos_parameters eos,CCTK_REAL rho_in, CCTK_REAL &P_cold,CCTK_REAL &eps_cold );
void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold( CCTK_REAL *U, const igm_eos_parameters eos,
                                                                   CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                   CCTK_REAL &Gamma_cold );
void compute_entropy_function( const igm_eos_parameters eos,
                               const CCTK_REAL rho,
                               const CCTK_REAL P,
                               CCTK_REAL *restrict S );

void reset_prims_to_atmosphere( const igm_eos_parameters eos, CCTK_REAL *restrict PRIMS );
//------------------------------------------------

//---------- Tabulated Equation of State ---------
void initialize_Tabulated_EOS_parameters_from_input( const CCTK_REAL cctk_time, igm_eos_parameters& eos );

void compute_remaining_prims_on_right_and_left_face( const igm_eos_parameters eos,
                                                     const cGH *restrict cctkGH,
                                                     const CCTK_INT *restrict cctk_lsh,
                                                     const gf_and_gz_struct *restrict in_prims,
                                                     gf_and_gz_struct *restrict out_prims_r,
                                                     gf_and_gz_struct *restrict out_prims_l );

void enforce_table_bounds_rho_Ye_T( const igm_eos_parameters& eos,
                                    CCTK_REAL *restrict rho,
                                    CCTK_REAL *restrict Ye,
                                    CCTK_REAL *restrict T );

void enforce_table_bounds_rho_Ye_eps( const igm_eos_parameters& eos,
                                      CCTK_REAL *restrict rho,
                                      CCTK_REAL *restrict Ye,
                                      CCTK_REAL *restrict eps );

void enforce_table_bounds_rho_Ye_S( const igm_eos_parameters& eos,
                                    CCTK_REAL *restrict rho,
                                    CCTK_REAL *restrict Ye,
                                    CCTK_REAL *restrict S );

void enforce_table_bounds_rho_Ye_P( const igm_eos_parameters& eos,
                                    CCTK_REAL *restrict rho,
                                    CCTK_REAL *restrict Ye,
                                    CCTK_REAL *restrict P );
//------------------------------------------------

#endif // IGM_EOS_HEADERS_HH


