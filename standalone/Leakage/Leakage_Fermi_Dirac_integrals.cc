// Returns Fermi-Dirac integral according to
// [1] Takahashi, Eid, and Hillebrandt (1978)
// https://adsabs.harvard.edu/pdf/1978A%26A....67..185T
#include "Leakage.hh"

REAL fermi_dirac_integral(const int N, const double z) {
  if( z > 1e-3 )
    switch(N) {
    case(0):
      // https://www.wolframalpha.com/input/?i=integral+of+1%2F%28exp%28x-z%29%2B1%29dx+from+0+to+infinity
      return(log(1.0+exp(z)));
      break;
    case(1):
      // First equation in Eqs. (A2) in [1].
      return((0.5*z*z+1.6449)/(1.0+exp(-1.6855*z)));
      break;
    case(2):
      // Second equation in Eqs. (A2) in [1].
      return((z*z*z/3.0 + 3.2899*z)/(1.0-exp(-1.8246*z)));
      break;
    case(3):
      // Third equation in Eqs. (A2) in [1].
      return((z*z*z*z/4.0 + 4.9348*z*z + 11.3644)/(1.0+exp(-1.9039*z)));
      break;
    case(4):
      // Fourth equation in Eqs. (A2) in [1].
      return((z*z*z*z*z/5.0 + 6.5797*z*z*z + 45.4576*z)/(1.0-exp(-1.9484*z)));
      break;
    case(5):
      // Fifth equation in Eqs. (A2) in [1].
      return((z*z*z*z*z*z/6.0 + 8.2247*z*z*z*z + 113.6439*z*z + 236.5323)/(1.0+exp(-1.9727*z)));
      break;
  }
  else {
    switch(N) {
    case(0):
      // https://www.wolframalpha.com/input/?i=integral+of+1%2F%28exp%28x-z%29%2B1%29dx+from+0+to+infinity
      return(log(1.0+exp(z)));
      break;
    case(1):
      // First equation in Eqs. (A3) in [1].
      return(exp(z) / (1.0 + 0.2159 * exp(0.8857 * z)));
      break;
    case(2):
      // Second equation in Eqs. (A3) in [1].
      return((2.0 * exp(z)) / (1.0 + 0.1092 * exp(0.8908 * z)));
      break;
    case(3):
      // Third equation in Eqs. (A3) in [1].
      return((6.0 * exp(z)) / (1.0 + 0.0559 * exp(0.9069 * z)));
      break;
    case(4):
      // Fourth equation in Eqs. (A3) in [1].
      return((24.0 * exp(z)) / (1.0 + 0.0287 * exp(0.9257 * z)));
      break;
    case(5):
      // Fifth equation in Eqs. (A3) in [1].
      return((120.0 * exp(z)) / (1.0 + 0.0147 * exp(0.9431 * z)));
      break;
    }
  }
}
