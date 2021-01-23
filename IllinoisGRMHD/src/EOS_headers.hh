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

  // What EOS will we use? IllinoisGRMHD
  // supports Hybrid or Tabulated (nuc_eos).
  CCTK_INT key;

  // Atmospheric density
  CCTK_REAL rho_atm;

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
  //------------------------------------------------

} igm_eos_parameters;
//------------------------------------------------

//----------- Hybrid Equation of State -----------
void print_EOS_Hybrid( igm_eos_parameters eos );
void setup_K_ppoly_tab__and__eps_integ_consts( igm_eos_parameters &eos );
void initialize_igm_eos_parameters_from_input( igm_eos_parameters &eos );
int find_polytropic_K_and_Gamma_index( igm_eos_parameters eos, CCTK_REAL rho_in );
int find_polytropic_K_and_Gamma_index_from_P( const igm_eos_parameters eos, const CCTK_REAL P_in );
void compute_P_cold__eps_cold(igm_eos_parameters eos,CCTK_REAL rho_in, CCTK_REAL &P_cold,CCTK_REAL &eps_cold );
void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold( CCTK_REAL *U, igm_eos_parameters &eos, CCTK_REAL Gamma_th,
                                                                   CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                   CCTK_REAL &Gamma_cold );
//------------------------------------------------

#endif // __EOS_HEADERS__


