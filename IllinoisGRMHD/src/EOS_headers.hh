// Thorn      : IllinoisGRMHD
// File       : EOS_headers.hh
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file contains function prototypes for the
//              EOS functions available in IllinoisGRMHD

#ifndef __EOS_HEADERS__
#define __EOS_HEADERS__

//----------- Hybrid Equation of State -----------
void print_EOS_Hybrid( eos_struct eos );
void setup_K_ppoly_tab__and__eps_integ_consts( eos_struct &eos );
void initialize_EOS_struct_from_input( eos_struct &eos );
int find_polytropic_K_and_Gamma_index( eos_struct eos, CCTK_REAL rho_in );
int find_polytropic_K_and_Gamma_index_from_P( const eos_struct eos, const CCTK_REAL P_in );
void compute_P_cold__eps_cold(eos_struct eos,CCTK_REAL rho_in, CCTK_REAL &P_cold,CCTK_REAL &eps_cold );
void compute_P_cold__eps_cold__dPcold_drho__eps_th__h__Gamma_cold( CCTK_REAL *U, eos_struct &eos, CCTK_REAL Gamma_th,
                                                                   CCTK_REAL &P_cold,CCTK_REAL &eps_cold,CCTK_REAL &dPcold_drho,CCTK_REAL &eps_th,CCTK_REAL &h,
                                                                   CCTK_REAL &Gamma_cold );
//------------------------------------------------

#endif // __EOS_HEADERS__


