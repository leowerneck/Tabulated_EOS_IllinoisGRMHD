// Thorn      : IllinoisGRMHD
// File       : EOS_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains function prototypes for the
//              EOS functions available in IllinoisGRMHD

#ifndef __EOS_HEADERS__
#define __EOS_HEADERS__

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
  // Which variable to reconstruct in PPM
  CCTK_INT PPM_reconstructed_var;
  // Whether or not to evolve the entropy
  bool evolve_entropy;
  // Atmospheric density
  CCTK_REAL rho_b_atm;
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
  // Number 
  CCTK_REAL eps_integ_const[MAX_EOS_PARAMS];
  //------------------------------------------------

  //---------- Tabulated Equation of State ---------
  // Atmospheric electron fraction
  CCTK_REAL Ye_atm;
  // Atmospheric temperature
  CCTK_REAL T_atm;
  // Atmospheric pressure
  CCTK_REAL P_atm;
  // Atmospheric specific internal energy
  CCTK_REAL eps_atm;
  // Atmospheric entropy
  CCTK_REAL S_atm;
  // Root-finding precision for table inversions
  CCTK_REAL root_finding_precision;
  //------------------------------------------------

} igm_eos_parameters;
//------------------------------------------------

void initialize_igm_eos_parameters_from_input( const int* igm_eos_key, igm_eos_parameters &eos );

//----------- Hybrid Equation of State -----------
void print_EOS_Hybrid( igm_eos_parameters eos );
void setup_K_ppoly_tab__and__eps_integ_consts( igm_eos_parameters &eos );
void initialize_Hybrid_EOS_parameters_from_input( igm_eos_parameters &eos );
int find_polytropic_K_and_Gamma_index( igm_eos_parameters eos, CCTK_REAL rho_in );
int find_polytropic_K_and_Gamma_index_from_P( const igm_eos_parameters eos, const CCTK_REAL P_in );
void compute_P_cold__eps_cold(igm_eos_parameters eos,CCTK_REAL rho_in, CCTK_REAL &P_cold,CCTK_REAL &eps_cold );
void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold( CCTK_REAL *U, igm_eos_parameters &eos, CCTK_REAL Gamma_th,
                                                                   CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                   CCTK_REAL &Gamma_cold );
//------------------------------------------------

//---------- Tabulated Equation of State ---------
void initialize_Tabulated_EOS_parameters_from_input( igm_eos_parameters& eos );

void get_P_and_eps_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                      const CCTK_REAL rho,
                                      const CCTK_REAL Y_e,
                                      const CCTK_REAL T,
                                      CCTK_REAL *restrict P,
                                      CCTK_REAL *restrict eps );

void get_P_eps_and_S_from_rho_Ye_and_T( const igm_eos_parameters eos,
                                        const CCTK_REAL rho,
                                        const CCTK_REAL Y_e,
                                        const CCTK_REAL T,
                                        CCTK_REAL *restrict P,
                                        CCTK_REAL *restrict eps,
                                        CCTK_REAL *restrict S );
//------------------------------------------------

#endif // __EOS_HEADERS__


