// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_validation.cc
// Author(s)  : Leo Werneck (wernecklr@gmail.com)
// Description: This file validates the implementation of this thorn
//              by comparing results from the functions it provides
//              against the ones provided by EOS_Omni.

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "WVU_EOS_Tabulated_headers.hh"

extern "C"
void WVU_EOS_validation_against_EOS_Omni(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS_WVU_EOS_validation_against_EOS_Omni;
  DECLARE_CCTK_PARAMETERS;
  
  // --------------------------------------------
  // --------- Basic information output ---------
  // --------------------------------------------
  CCTK_VInfo(CCTK_THORNSTRING,".----------------------------------------.");
  CCTK_VInfo(CCTK_THORNSTRING,"| Unit test: validation against EOS_Omni |");
  CCTK_VInfo(CCTK_THORNSTRING,".----------------------------------------.");
  CCTK_VInfo(CCTK_THORNSTRING,"Failure of any of the tests below can lead to major problems, thus we will terminate the run instead.\n");
  
  int max_length = 0;
  for(int i=0;i<WVU_EOS::ntables;i++) {
    int aux_length = WVU_EOS::table_var_names[i].length();
    if( aux_length > max_length ) max_length = aux_length;
  }

  for(int i=0;i<WVU_EOS::ntables;i++) {
    while( (int)WVU_EOS::table_var_names[i].length() < max_length ) {
      WVU_EOS::table_var_names[i] += " ";
    }
    CCTK_VInfo(CCTK_THORNSTRING,"The EOS table key for variable %s is %d",WVU_EOS::table_var_names[i].c_str(),i);
  }
  printf("\n");
  // --------------------------------------------

  // Set values of rho, Ye, and T for unit tests
  CCTK_INT unit_test_npoints = 4;

  // Minimum/maximum values of rho, Ye, T that will be tested
  CCTK_REAL unit_test_lr_max = log(0.9*nuc_eos::eos_rhomax);
  CCTK_REAL unit_test_lr_min = log(1.1*nuc_eos::eos_rhomin);

  CCTK_REAL unit_test_lt_max = log(0.9*nuc_eos::eos_tempmax);
  CCTK_REAL unit_test_lt_min = log(1.1*nuc_eos::eos_tempmin);

  CCTK_REAL unit_test_ly_max = log(0.9*nuc_eos::eos_yemax);
  CCTK_REAL unit_test_ly_min = log(1.1*nuc_eos::eos_yemin);

  // Compute step sizes
  CCTK_REAL dlr = (unit_test_lr_max - unit_test_lr_min)/unit_test_npoints;
  CCTK_REAL dlt = (unit_test_lt_max - unit_test_lt_min)/unit_test_npoints;
  CCTK_REAL dly = (unit_test_ly_max - unit_test_ly_min)/unit_test_npoints;

  // Basic test quantities
  const CCTK_INT eoskey   = EOS_Omni_GetHandle("nuc_eos");
  const CCTK_INT npoints  = 1;
  const CCTK_INT haveeps  = 0;
  const CCTK_INT havetemp = 1;
  const CCTK_INT haveent  = 2;
  const CCTK_REAL prec    = 1e-10;

  // -------------------------------------------------------------------
  // ----- Unit test for function: WVU_EOS_P_and_eps_from_rho_Ye_T -----
  // -------------------------------------------------------------------
  {
    CCTK_INT counter = 0;
    CCTK_INT passed1 = 0;
    CCTK_INT passed2 = 0;
    CCTK_INT passed3 = 0;
    CCTK_INT passed4 = 0;
    for(int k=0;k<=unit_test_npoints;k++) {
      for(int j=0;j<=unit_test_npoints;j++) {
        for(int i=0;i<=unit_test_npoints;i++) {
          // Set rho, Ye, and T
          CCTK_REAL xrho  = exp(unit_test_lr_min + i*dlr);
          CCTK_REAL xtemp = exp(unit_test_lt_min + j*dlt);
          CCTK_REAL xye   = exp(unit_test_ly_min + k*dly);
          // Compute P and eps using WVU_EOS
          CCTK_REAL xprs_wvu  = 0.0;
          CCTK_REAL xeps_wvu  = 0.0;
          WVU_EOS_P_and_eps_from_rho_Ye_T( xrho,xye,xtemp, &xprs_wvu,&xeps_wvu );
          // Compute P and eps using EOS_Omni
          CCTK_INT anyerr     = 0;
          CCTK_INT keyerr     = 0;
          CCTK_REAL xprs_omni = 0.0;
          CCTK_REAL xeps_omni = 0.0;
          EOS_Omni_press( eoskey,havetemp,prec,npoints,
                          &xrho,&xeps_omni,&xtemp,&xye,&xprs_omni,
                          &keyerr,&anyerr );

          // Compute the relative error between the results
          CCTK_REAL prs_err = fabs(1.0 - xprs_wvu/xprs_omni);
          CCTK_REAL eps_err = fabs(1.0 - xeps_wvu/xeps_omni);
          // Check if results are correct
          if( (prs_err < 1e-13) && (eps_err < 1e-13) ) {
            passed1++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_and_eps_from_rho_Ye_T FAILED! %e %e %e -> WVU: %e %e | Omni: %e %e. Relative errors: %e %e. ABORTING.",
                        xrho,xye,xtemp,xprs_wvu,xeps_wvu,xprs_omni,xeps_omni,prs_err,eps_err);
          }
          // Now that we already have P and eps, let us test other functions
          // Compute P and T using WVU_EOS
          xprs_wvu             = 0.0;
          CCTK_REAL xtemp_wvu  = xtemp;
          WVU_EOS_P_and_T_from_rho_Ye_eps( xrho,xye,xeps_wvu, &xprs_wvu,&xtemp_wvu );
          // Compute P and T using EOS_Omni
          xprs_omni            = 0.0;
          CCTK_REAL xtemp_omni = xtemp;
          CCTK_REAL xent_omni  = 0.0;
          CCTK_REAL xdummy     = 0.0;
          EOS_Omni_short( eoskey,haveeps,prec,npoints,
                          &xrho,&xeps_omni,&xtemp_omni,&xye,&xprs_omni,&xent_omni,
                          &xdummy,&xdummy,&xdummy,&xdummy,&xdummy,
                          &keyerr,&anyerr );

          // Compute the relative error between the results
          prs_err            = fabs(1.0 -  xprs_wvu/xprs_omni);
          CCTK_REAL temp_err = fabs(1.0 - xtemp_wvu/xtemp_omni);
          // Check if results are correct
          if( (prs_err < 1e-13) && (temp_err < 1e-13) ) {
            passed2++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_and_T_from_rho_Ye_eps FAILED! %e %e %e -> WVU: %e %e | Omni: %e %e. Relative errors: %e %e. ABORTING.",
                        xrho,xye,xeps_wvu,xprs_wvu,xtemp_wvu,xprs_omni,xtemp_omni,prs_err,temp_err);
          }

          // Compute P, S, and T(eps) using WVU_EOS
          xprs_wvu           = 0.0;
          xtemp_wvu          = xtemp;
          CCTK_REAL xent_wvu = 0.0;
          WVU_EOS_P_S_and_T_from_rho_Ye_eps( xrho,xye,xeps_wvu,&xprs_wvu,&xent_wvu,&xtemp_wvu );
          // Now compute the errors
                    temp_err          = fabs(1.0 - xtemp_wvu/xtemp_omni);
          CCTK_REAL ent_err = fabs(1.0 - xent_wvu /xent_omni);
          if( (prs_err < 1e-13) && (temp_err < 1e-13) && (ent_err < 1e-13) ) {
            passed3++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_S_and_T_from_rho_Ye_eps FAILED! %e %e %e -> WVU: %e %e %e | Omni: %e %e %e. Relative errors: %e %e %e. ABORTING.",
                        xrho,xye,xeps_wvu,xprs_wvu,xent_wvu,xtemp_wvu,xprs_omni,xent_omni,xtemp_omni,prs_err,ent_err,temp_err);
          }

          // Compute P, eps, T(S) using WVU_EOS
          xprs_wvu  = 0.0;
          xeps_wvu  = 0.0;
          xtemp_wvu = xtemp;
          WVU_EOS_P_eps_and_T_from_rho_Ye_S( xrho,xye,xent_wvu, &xprs_wvu,&xeps_wvu,&xtemp_wvu );
          // Compute P, eps, T(S) using EOS_Omni
          EOS_Omni_short( eoskey,haveent,prec,npoints,
                          &xrho,&xeps_omni,&xtemp_omni,&xye,&xprs_omni,&xent_omni,
                          &xdummy,&xdummy,&xdummy,&xdummy,&xdummy,
                          &keyerr,&anyerr );

          // Now compute the errors
          prs_err  = fabs(1.0 - xprs_wvu /xprs_omni);
          eps_err  = fabs(1.0 - xeps_wvu /xeps_omni);
          temp_err = fabs(1.0 - xtemp_wvu/xtemp_omni);
          if( (prs_err < 1e-13) && (eps_err < 1e-13) && (temp_err < 1e-13) ) {
            passed4++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_eps_and_T_from_rho_Ye_S FAILED! %e %e %e -> WVU: %e %e %e | Omni: %e %e %e. Relative errors: %e %e %e. ABORTING.",
                        xrho,xye,xent_wvu,xprs_wvu,xeps_wvu,xtemp_wvu,xprs_omni,xeps_omni,xtemp_omni,prs_err,eps_err,temp_err);
          }
          counter++;
        }
      }
    }
    // If we got here then all tests above passed
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_and_eps_from_rho_Ye_T                            : %d/%d tests PASSED",passed1,counter);
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_and_T_from_rho_Ye_eps                            : %d/%d tests PASSED",passed2,counter);
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_S_and_T_from_rho_Ye_eps                          : %d/%d tests PASSED",passed3,counter);
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_eps_and_T_from_rho_Ye_S                          : %d/%d tests PASSED",passed4,counter);
  }

  // ---------------------------------------------------------------------
  // ----- Unit test for function: WVU_EOS_P_eps_and_S_from_rho_Ye_T -----
  // ---------------------------------------------------------------------
  {
    CCTK_INT counter = 0;
    CCTK_INT passed  = 0;
    for(int k=0;k<=unit_test_npoints;k++) {
      for(int j=0;j<=unit_test_npoints;j++) {
        for(int i=0;i<=unit_test_npoints;i++) {
          // Set rho, Ye, and T
          CCTK_REAL xrho  = exp(unit_test_lr_min + i*dlr);
          CCTK_REAL xtemp = exp(unit_test_lt_min + j*dlt);
          CCTK_REAL xye   = exp(unit_test_ly_min + k*dly);
          // Compute P and eps using WVU_EOS
          CCTK_REAL xprs_wvu  = 0.0;
          CCTK_REAL xeps_wvu  = 0.0;
          CCTK_REAL xent_wvu  = 0.0;
          WVU_EOS_P_eps_and_S_from_rho_Ye_T( xrho,xye,xtemp, &xprs_wvu,&xeps_wvu,&xent_wvu );
          // Compute P and eps using EOS_Omni
          CCTK_INT anyerr     = 0;
          CCTK_INT keyerr     = 0;
          CCTK_REAL xprs_omni = 0.0;
          CCTK_REAL xeps_omni = 0.0;
          CCTK_REAL xent_omni = 0.0;
          CCTK_REAL xdummy    = 0.0;
          EOS_Omni_short( eoskey,havetemp,prec,npoints,
                          &xrho,&xeps_omni,&xtemp,&xye,&xprs_omni,&xent_omni,
                          &xdummy,&xdummy,&xdummy,&xdummy,&xdummy,
                          &keyerr,&anyerr );
          // Compute the relative error between the results
          CCTK_REAL prs_err = fabs(1.0 - xprs_wvu/xprs_omni);
          CCTK_REAL eps_err = fabs(1.0 - xeps_wvu/xeps_omni);
          CCTK_REAL ent_err = fabs(1.0 - xent_wvu/xent_omni);
          // Check if results are correct
          if( (prs_err < 1e-13) && (eps_err < 1e-13) && (ent_err < 1e-13) ) {
            passed++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_eps_and_S_from_rho_Ye_T FAILED! %e %e %e -> WVU: %e %e %e | Omni: %e %e %e. Relative errors: %e %e %e. ABORTING.",
                        xrho,xye,xtemp,xprs_wvu,xeps_wvu,xent_wvu,xprs_omni,xeps_omni,xent_omni,prs_err,eps_err,ent_err);
          }
          counter++;
        }
      }
    }
    // If we got here then all tests above passed
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_eps_and_S_from_rho_Ye_T                          : %d/%d tests PASSED",counter,counter);
  }

  // ---------------------------------------------------------------------
  // --- Unit test for function: WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T ---
  // ---------------------------------------------------------------------
  {
    CCTK_INT counter = 0;
    CCTK_INT passed  = 0;
    for(int k=0;k<=unit_test_npoints;k++) {
      for(int j=0;j<=unit_test_npoints;j++) {
        for(int i=0;i<=unit_test_npoints;i++) {
          // Set rho, Ye, and T
          CCTK_REAL xrho  = exp(unit_test_lr_min + i*dlr);
          CCTK_REAL xtemp = exp(unit_test_lt_min + j*dlt);
          CCTK_REAL xye   = exp(unit_test_ly_min + k*dly);
          // Compute P and eps using WVU_EOS
          CCTK_REAL xprs_wvu = 0.0;
          CCTK_REAL xeps_wvu = 0.0;
          CCTK_REAL xent_wvu = 0.0;
          CCTK_REAL xcs2_wvu = 0.0;
          WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T( xrho,xye,xtemp, &xprs_wvu,&xeps_wvu,&xent_wvu,&xcs2_wvu );
          // Compute P and eps using EOS_Omni
          CCTK_INT anyerr     = 0;
          CCTK_INT keyerr     = 0;
          CCTK_REAL xprs_omni = 0.0;
          CCTK_REAL xeps_omni = 0.0;
          CCTK_REAL xent_omni = 0.0;
          CCTK_REAL xcs2_omni = 0.0;
          CCTK_REAL xdummy    = 0.0;
          EOS_Omni_short( eoskey,havetemp,prec,npoints,
                          &xrho,&xeps_omni,&xtemp,&xye,&xprs_omni,&xent_omni,
                          &xcs2_omni,&xdummy,&xdummy,&xdummy,&xdummy,
                          &keyerr,&anyerr );
          // Compute the relative error between the results
          CCTK_REAL prs_err = fabs(1.0 - xprs_wvu/xprs_omni);
          CCTK_REAL eps_err = fabs(1.0 - xeps_wvu/xeps_omni);
          CCTK_REAL ent_err = fabs(1.0 - xent_wvu/xent_omni);
          CCTK_REAL cs2_err = fabs(1.0 - xcs2_wvu/xcs2_omni);
          // Check if results are correct
          if( (prs_err < 1e-13) && (eps_err < 1e-13) && (ent_err < 1e-13) && (cs2_err < 1e-13) ) {
            passed++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T FAILED! %e %e %e -> WVU: %e %e %e %e | Omni: %e %e %e %e. Relative errors: %e %e %e %e. ABORTING.",
                        xrho,xye,xtemp,xprs_wvu,xeps_wvu,xent_wvu,xcs2_wvu,xprs_omni,xeps_omni,xent_omni,xcs2_omni,prs_err,eps_err,ent_err,cs2_err);
          }
          counter++;
        }
      }
    }
    // If we got here then all tests above passed
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T                      : %d/%d tests PASSED",counter,counter);
  }

  // ----------------------------------------------------------------------
  // --- Unit test for function: WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T ---
  // ----------------------------------------------------------------------
  {
    CCTK_INT counter = 0;
    CCTK_INT passed  = 0;
    for(int k=0;k<=unit_test_npoints;k++) {
      for(int j=0;j<=unit_test_npoints;j++) {
        for(int i=0;i<=unit_test_npoints;i++) {
          // Set rho, Ye, and T
          CCTK_REAL xrho  = exp(unit_test_lr_min + i*dlr);
          CCTK_REAL xtemp = exp(unit_test_lt_min + j*dlt);
          CCTK_REAL xye   = exp(unit_test_ly_min + k*dly);
          // Compute P and eps using WVU_EOS
          CCTK_REAL xprs_wvu  = 0.0;
          CCTK_REAL xeps_wvu  = 0.0;
          CCTK_REAL xdedt_wvu = 0.0;
          WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T( xrho,xye,xtemp, &xprs_wvu,&xeps_wvu,&xdedt_wvu );
          // Compute P and eps using EOS_Omni
          CCTK_INT anyerr      = 0;
          CCTK_INT keyerr      = 0;
          CCTK_REAL xprs_omni  = 0.0;
          CCTK_REAL xeps_omni  = 0.0;
          CCTK_REAL xdedt_omni = 0.0;
          CCTK_REAL xdummy     = 0.0;
          EOS_Omni_short( eoskey,havetemp,prec,npoints,
                          &xrho,&xeps_omni,&xtemp,&xye,&xprs_omni,&xdummy,
                          &xdummy,&xdedt_omni,&xdummy,&xdummy,&xdummy,
                          &keyerr,&anyerr );
          // Compute the relative error between the results
          CCTK_REAL prs_err  = fabs(1.0 -  xprs_wvu/xprs_omni);
          CCTK_REAL eps_err  = fabs(1.0 -  xeps_wvu/xeps_omni);
          CCTK_REAL dedt_err = fabs(1.0 - xdedt_wvu/xdedt_omni);
          // Check if results are correct
          if( (prs_err < 1e-13) && (eps_err < 1e-13) && (dedt_err < 1e-13) ) {
            passed++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_eps_S_and_cs2_from_rho_Ye_T FAILED! %e %e %e -> WVU: %e %e %e | Omni: %e %e %e. Relative errors: %e %e %e. ABORTING.",
                        xrho,xye,xtemp,xprs_wvu,xeps_wvu,xdedt_wvu,xprs_omni,xeps_omni,xdedt_omni,prs_err,eps_err,dedt_err);
          }
          counter++;
        }
      }
    }
    // If we got here then all tests above passed
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_eps_and_depsdT_from_rho_Ye_T                     : %d/%d tests PASSED",counter,counter);
  }

  // -------------------------------------------------------------------------------------------
  // --- Unit test for function: WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T ---
  // -------------------------------------------------------------------------------------------
  {
    CCTK_INT counter = 0;
    CCTK_INT passed  = 0;
    for(int k=0;k<=unit_test_npoints;k++) {
      for(int j=0;j<=unit_test_npoints;j++) {
        for(int i=0;i<=unit_test_npoints;i++) {
          // Set rho, Ye, and T
          CCTK_REAL xrho  = exp(unit_test_lr_min + i*dlr);
          CCTK_REAL xtemp = exp(unit_test_lt_min + j*dlt);
          CCTK_REAL xye   = exp(unit_test_ly_min + k*dly);
          // Compute P and eps using WVU_EOS
          CCTK_REAL xprs_wvu  = 0.0;
          CCTK_REAL xeps_wvu  = 0.0;
          CCTK_REAL xdPdr_wvu = 0.0;
          CCTK_REAL xdPdt_wvu = 0.0;
          CCTK_REAL xdedr_wvu = 0.0;
          CCTK_REAL xdedt_wvu = 0.0;
          WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T( xrho,xye,xtemp, &xprs_wvu,&xeps_wvu,
                                                                       &xdPdr_wvu,&xdPdt_wvu,&xdedr_wvu,&xdedt_wvu );
          // Compute P and eps using EOS_Omni
          CCTK_INT anyerr      = 0;
          CCTK_INT keyerr      = 0;
          CCTK_REAL xprs_omni  = 0.0;
          CCTK_REAL xeps_omni  = 0.0;
          CCTK_REAL xdPdr_omni = 0.0;
          CCTK_REAL xdPde_omni = 0.0;
          CCTK_REAL xdPdt_omni = 0.0;
          CCTK_REAL xdedr_omni = 0.0;
          CCTK_REAL xdedt_omni = 0.0;
          CCTK_REAL xdummy     = 0.0;
          EOS_Omni_short( eoskey,havetemp,prec,npoints,
                          &xrho,&xeps_omni,&xtemp,&xye,&xprs_omni,&xdummy,
                          &xdummy,&xdedt_omni,&xdPde_omni,&xdPdr_omni,&xdummy,
                          &keyerr,&anyerr );
          xdPdt_omni = (xdPde_omni)*(xdedt_omni);
          xdedr_omni = (xdPdr_omni)/(xdPde_omni);
          // Compute the relative error between the results
          CCTK_REAL prs_err  = fabs(1.0 -  xprs_wvu/xprs_omni);
          CCTK_REAL eps_err  = fabs(1.0 -  xeps_wvu/xeps_omni);
          CCTK_REAL dPdr_err = fabs(1.0 - xdPdr_wvu/xdPdr_omni);
          CCTK_REAL dPdt_err = fabs(1.0 - xdPdt_wvu/xdPdt_omni);
          CCTK_REAL dedr_err = fabs(1.0 - xdedr_wvu/xdedr_omni);
          CCTK_REAL dedt_err = fabs(1.0 - xdedt_wvu/xdedt_omni);
          // Check if results are correct
          if( (prs_err  < 1e-13) && (eps_err  < 1e-13) && (dedt_err < 1e-13) &&
              (dPdr_err < 1e-13) && (dPdt_err < 1e-13) && (dedr_err < 1e-13) ) {
            passed++;
          }
          else {
            CCTK_VError(__LINE__,__FILE__,CCTK_THORNSTRING,
                        "Unit test for function WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T FAILED! %e %e %e -> WVU: %e %e %e %e %e %e | Omni: %e %e %e %e %e %e. Relative errors: %e %e %e %e %e %e. ABORTING.",
                        xrho,xye,xtemp,xprs_wvu,xeps_wvu,xdPdr_wvu,xdPdt_wvu,xdedr_wvu,xdedt_wvu,
                        xprs_omni,xeps_omni,xdPdr_wvu,xdPdt_wvu,xdedr_wvu,xdedt_wvu,
                        prs_err,eps_err,dPdr_err,dPdt_err,dedr_err,dedt_err);
          }
          counter++;
        }
      }
    }
    // If we got here then all tests above passed
    CCTK_VInfo(CCTK_THORNSTRING,"Unit test: WVU_EOS_P_eps_dPdrho_dPdT_depsdrho_and_depsdT_from_rho_Ye_T: %d/%d tests PASSED",counter,counter);
  }

  CCTK_VInfo(CCTK_THORNSTRING,"All tests PASSED!");
  exit(1);
  
}
