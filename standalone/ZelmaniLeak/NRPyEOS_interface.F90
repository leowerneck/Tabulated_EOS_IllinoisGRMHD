subroutine NRPyEOS_fortran_interface(eos_params, rho, Y_e, T, P, eps)

  use, intrinsic :: iso_c_binding, only: c_double
  use NRPyEOS

  implicit none

  interface
     subroutine NRPyEOS_P_and_eps_from_rho_Ye_T_fortran(eos_params, rho, Y_e, T, P, eps) BIND(C)
       use, intrinsic :: iso_c_binding, only: c_double
       use NRPyEOS
       implicit none
       type(NRPyEOS_params), intent(in) :: eos_params
       real(kind=c_double), intent(in)  :: rho, Y_e, T
       real(kind=c_double), intent(out) :: P, eps
     end subroutine NRPyEOS_P_and_eps_from_rho_Ye_T_fortran
  end interface

  type(NRPyEOS_params), intent(in) :: eos_params
  real(kind=c_double), intent(in) :: rho, Y_e, T
  real(kind=c_double), intent(out) :: P, eps

  call NRPyEOS_P_and_eps_from_rho_Ye_T_fortran(eos_params, rho, Y_e, T, P, eps)

end subroutine NRPyEOS_fortran_interface
