#ifndef LEAKAGE_TESTS_HH_
#define LEAKAGE_TESTS_HH_

//******************************************************************************
// These functions can be found in the files with the same
// name inside the optically_thin subdirectory
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
//******************************************************************************

//******************************************************************************
// These functions can be found in the files with the same
// name inside the constant_density_sphere subdirectory
void ConstantDensitySphere_NRPyLeakage();

void ConstantDensitySphere_harm_leakage();
//******************************************************************************

#endif // LEAKAGE_TESTS_HH_
