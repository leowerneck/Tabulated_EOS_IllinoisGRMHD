#ifndef LEAKAGE_TESTS_HH_
#define LEAKAGE_TESTS_HH_

void OpticallyThinGas_NRPyLeakage( const int which_constants,
                                   const NRPyEOS_params *restrict eos_params,
                                   const REAL initial_rho_b,
                                   const REAL initial_Y_e,
                                   const REAL initial_T );

void OpticallyThinGas_harm_leakage( const int which_constants,
                                    const NRPyEOS_params *restrict eos_params,
                                    const REAL initial_rho_b,
                                    const REAL initial_Y_e,
                                    const REAL initial_T );

#endif // LEAKAGE_TESTS_HH_
