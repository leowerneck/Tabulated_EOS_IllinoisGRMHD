#ifndef NEUTRINOS
#define NEUTRINOS
//#include "decs.h"
#include "NRPyEOS.h"

// HARM macros (from macros.h)
#define NDIM     (4) /* number of spacetime dimensions, 0=time */
#define RHO      (0) /* rest-mass density                      */
#define UU       (1) /* internal energy density                */ 
#define U1       (2) /* spatial velocity component             */
#define U2       (3) /* spatial velocity component             */
#define U3       (4) /* spatial velocity component             */
#define B1       (5) /* spatial magnetic field component       */
#define B2       (6) /* spatial magnetic field component       */
#define B3       (7) /* spatial magnetic field component       */
#define YE       (8) /* electron fraction component            */
#define TEMP     (9) /* Temperature (in MeV)                   */
#define NUMPRIMS (10)
#define N_OPTICAL_DEPTHS (6)

// HARM types (from globals.h)
struct of_geom {
  double gcon[NDIM][NDIM];
  double gcov[NDIM][NDIM];
  double g;
  double g_inv; 
  double alpha;
  double beta[NDIM-1];
  double ncon[NDIM];
};

struct of_state {
  double ucon[NDIM];
  double ucov[NDIM];
  double bcon[NDIM];
  double bcov[NDIM];
  double bsq;
  double p; 
  double eta;
  double ptot;
};

/*
 * Particle types
 */
#define nu__e       (0)
#define nu__e_bar   (1)
#define nu__mu      (2)
#define nu__mu_bar  (3)
#define nu__tau     (4)
#define nu__tau_bar (5)

#define NEUTRON     (6)
#define PROTON      (7)
#define ELECTRON    (8)
#define POSITRON    (9)
#define MUHAT       (10) 

#define NUMBER_OF_NEUTRINO_TYPES (6)

/*
 * Used in determining number fractions for electrons, protons, and neutrons
 */

#define COMPLETELY_DISSOCIATED          (0)

#define MATTER_DISTRIBUTION_ASSUMPTION   (COMPLETELY_DISSOCIATED)

/*
 * Used in deciding to compute neutrino-number transport or neutrino-energy transport
 */
#define NEUTRINO_NUMBER     (0)
#define NEUTRINO_ENERGY     (1)

/*
 * Used in switching between phase space blocking expressions
 */
#define SCATTERING                              (0)
#define BETA_PROCESS                            (1)
#define ELECTRON_POSITRON_PAIR_ANNIHILATION     (2)
#define PLASMON_DECAY                           (3)
#define NEUTRINO                                (4)
#define ANTI_NEUTRINO                           (5)

/*
 * Used to compute different energy moments
 */
#define ORDINARY    (0)
#define TILDE       (1)
#define STAR        (2)

/*
 * Used to define the Equation of state for the determination of Temperature
 */
#ifndef GAMMA
#define GAMMA   (0)
#endif

#define EOS     (GAMMA)

/*
 * ---- Constants ----
 *
 * h :=  6.6260693e-27
 *       Planck's constant, measured in ((cm^2 * g)/s)
 * c :=  2.99792458e10
 *       Speed of light, measured in (cm/s)
 * C_hc  :=  (h * c)
 * C_hc3 :=  (h * c)^3
 *
 * me := 9.1093826e-28
 *    Electron mass, measured in (g)
 * C_me2_c3  := me^2 * c^3
 * sigma_0   := 1.76e-44    after [1](A1) 
 * C_sigma_0__me2_c3 := sigma_0 / (me^2 * c^3)
 * A         := 6.02214086e23 1/mol
 C_A_sigma_0 =  A * sigma_0  
 C_alpha_neutrino_sq = (1.25)^2    after [1](A1) 
*/
#define C_hc    (1.9864456023253394e-16) //checked
#define C_hc3   (7.83844706784507644502e-48) //checked
#define C_me2_c3    (2.23583625945590876861e-23) //checked
#define C_sigma_0__me2_c3   (7.87177501285e-22) //ari: changed it, original: 6.4405431923285129459e-22
#define C_A_sigma_0  (1.05989679136e-20) // checked
#define C_alpha_neutrino_sq  (1.5625) //checked
#define C_fine_struct (0.00729735252051) //fine structure constant
#define C_me_c2_ev (0.5109989461e6) //me*c^2 in ev
#define C_hc3_ev (1.90589514992e-12) // where h is in eV*s
#define C_hc_ev (0.000123984197386)// where h is in eV*s
#define C_me2_c3_ev (8.71002308255) //in eV^2/(cm/s)
#define C_sigma_0__me2_c3_ev (2.02066054627e-45) //in cm^3/s 1/eV^2
#define C_sigma_0_cl (5.2763472608e-34) //sigma0*c
#define C_me2_c4 (2.61119922915e+11)


double get_Temperature(double energy_density, double rho, double electron_fraction);
double get_particle_fraction(int particle_type, double ye);
double get_degeneracy(const NRPyEOS_params *restrict eos_params,
                      int particle_type, double rho, double energy_density, 
                      double temperature, double electron_fraction);

void get_degeneracy_all(const NRPyEOS_params *restrict eos_params,
                        double rho, double energy_density, 
                        double temperature, double electron_fraction, 
                        double * eta_e, double * eta_p, double * eta_n, double * eta_muhat, 
                        double * eta__nu_e, double * eta__nu_e_bar,
                        double * eta__nu_mu, double * eta__nu_mu_bar, 
                        double * eta__nu_tau, double * eta__nu_tau_bar);

double get_Fermi_integral(int order, double eta);
double get_neutrino_number_density(int neutrino_flavor, double degeneracy, double temperature);
double get_neutrino_energy_density(int neutrino_flavor, double degeneracy, double temperature);
double get_energy_moment(int positron_or_electron, int symbol, double temperature, double electron_degeneracy);

double get_scattering_coefficient(int nucleon);
double get_phase_space_blocking_value(int interaction, int particle_type, double electron_degeneracy, double neutrino_degeneracy);
double get_spectrally_averaged_absorption_opacity(int neutrino_flavor, int nucleon, int type_of_transport,
                                                  double rho, double temperature, double electron_fraction,
                                                  double muhat_degeneracy, 
                                                  double electron_degeneracy, double neutrino_degeneracy);
double get_spectrally_averaged_scattering_opacity(int neutrino_flavor, int nucleon, int type_of_transport,
                                                  double rho, double temperature, double nucleon_fraction,
                                                  double nucleon_degeneracy, double neutrino_degeneracy);
double get_total_transport_opacity(const NRPyEOS_params *restrict eos_params,
                                   int neutrino_flavor, int type_of_transport, double rho, double energy_density, 
                                   double electron_fraction, double Temperature);

double beta_process_emission_rate__neutrino_number(int neutrino_flavor, double rho, double temperature,
                                                   double electron_fraction, double neutrino_degeneracy,
                                                   double muhat_degeneracy, double electron_degeneracy);
double beta_process_emission_rate__neutrino_energy(int neutrino_flavor, double rho, double temperature,
                                                   double electron_fraction, double neutrino_degeneracy,
                                                   double muhat_degeneracy, double electron_degeneracy);

double electron_positron_pair_annihilation_emission_rate__neutrino_number(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                                          double neutrino_degeneracy, double anti_neutrino_degeneracy);
double electron_positron_pair_annihilation_emission_rate__neutrino_energy(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                                          double neutrino_degeneracy, double anti_neutrino_degeneracy);

double plasmon_decay_emission_rate__neutrino_number(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                    double neutrino_degeneracy, double anti_neutrino_degeneracy);
double plasmon_decay_emission_rate__neutrino_energy(int neutrino_flavor, double temperature, double electron_degeneracy,
                                                    double neutrino_degeneracy, double anti_neutrino_degeneracy);



double get_optical_depth(int neutrino_flavor, int type_of_transport, double previous_optical_depth, double ds,
                         double neutrino_degeneracy, double rho, double energy_density, double electron_fraction,
                         double neutron_fraction, double proton_fraction, double electron_degeneracy,
                         double neutron_degeneracy, double proton_degeneracy);
double get_diffusion_timescale(const NRPyEOS_params *restrict eos_params,
                               int neutrino_flavor, int type_of_transport, double rho, 
                               double energy_density, double electron_fraction, double temperature,double Tau);
double get_inverse_emission_timescale(int neutrino_flavor, int type_of_transport, double rho, double temperature,
                                      double electron_fraction, double neutrino_degeneracy, double anti_neutrino_degeneracy,
                                      double muhat_degeneracy,
                                      double neutron_degeneracy, double proton_degeneracy, double electron_degeneracy);



double effective_emission_rate(const NRPyEOS_params *restrict eos_params,
                               int neutrino_flavor, int type_of_transport, 
                               double rho, double energy_density, double temperature,
                               double electron_fraction, double muhat_degeneracy,double neutron_fraction, double proton_fraction,
                               double neutrino_degeneracy, double anti_neutrino_degeneracy,
                               double neutron_degeneracy, double proton_degeneracy, double electron_degeneracy, double Tau);


double average_neutrino_energy_emitted(const NRPyEOS_params *restrict eos_params,
                                       int neutrino_flavor, double rho, 
                                       double energy_density, double temperature,
                                       double electron_fraction, double neutron_fraction, double proton_fraction,
                                       double Tau_number, double Tau_energy);

void neutrino_source_func(const NRPyEOS_params *restrict eos_params,
                          double *ph, double *optical_depth, struct of_state *q, struct of_geom *geom,
                          int ii, int jj, int kk, double *dU);

void neutrino_absorption_heating_rate(const NRPyEOS_params *restrict eos_params, double *ph, double *optical_depth,double *R_code_units, double *Q_code_units);

void neutrino_electron_Q(const NRPyEOS_params *restrict eos_params,
                         double *ph, double *optical_depth, double *Q_code_units );

void antineutrino_electron_Q(const NRPyEOS_params *restrict eos_params,
                             double *ph, double *optical_depth, double *Q_code_units );

void neutrino_x_Q(const NRPyEOS_params *restrict eos_params, double *ph, double *optical_depth, double *Q_code_units );

//typedef struct Nodes {
//    int origin;
//    int fence;

//    struct Nodes *next;

//} Node;

//typedef struct Node_List {
//    int size;
//    Node root;
//} List;

#endif      /* --- NEUTRINO_COOLING --- */
