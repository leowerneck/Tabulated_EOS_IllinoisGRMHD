#ifndef LEAKAGE_HH_
#define LEAKAGE_HH_

#include "Basic_defines.hh"

// Constants
namespace Leakage {
  // Defined below Eq. (A1) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  constexpr REAL sigma_0 = 1.76e-44;
  constexpr REAL sinthw2 = 0.23;
  constexpr REAL C_V     = 0.5 + 2*sinthw2;
  constexpr REAL C_A     = 0.5;
  constexpr REAL alpha   = 1.25;

  // Defined below Eq. (B9) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  constexpr REAL CV_minus_CA         = C_V - C_A;
  constexpr REAL CV_plus_CA          = C_V + C_A;
  constexpr REAL CV_minus_CA_sqrd    = CV_minus_CA*CV_minus_CA;
  constexpr REAL C1_plus_C2_nue_anue = CV_minus_CA_sqrd + CV_plus_CA*CV_plus_CA;

  // Defined below Eq. (B10) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  constexpr REAL CV_plus_CA_minus_2  = CV_plus_CA - 2.0;
  constexpr REAL C1_plus_C2_nux_anux = CV_minus_CA_sqrd + CV_plus_CA_minus_2*CV_plus_CA_minus_2;

  // Defined below Eq. (B12) in https://adsabs.harvard.edu/pdf/1996A%26A...311..532R
  constexpr REAL gamma_0 = 5.565e-2;

  // Physical constants
  constexpr REAL hc3 = 1.905895198207217e-30; // (h.c)^3 in (MeV.cm)^3
  constexpr REAL amu = 1.660539066600000e-24; // atomic mass unit in g
  constexpr REAL N_A = 6.022140760000000e+23; // Avogadro's number in 1/mol

  // Powers of Pi
  constexpr REAL pisqr = M_PI*M_PI;
}

// Leakage scheme function prototypes
REAL fermi_dirac_integral(const int N, const double z);

#endif // LEAKAGE_HH_
