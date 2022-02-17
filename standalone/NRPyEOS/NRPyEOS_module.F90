module NRPyEOS
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr
  type, bind(C) :: NRPyEOS_params
     integer(kind=c_int) :: nrho
     integer(kind=c_int) :: ntemp
     integer(kind=c_int) :: nye

     type(c_ptr) :: alltables
     type(c_ptr) :: epstable
     type(c_ptr) :: logrho
     type(c_ptr) :: logtemp
     type(c_ptr) :: yes

     real(kind=c_double) :: eos_rhomax , eos_rhomin
     real(kind=c_double) :: eos_tempmin, eos_tempmax
     real(kind=c_double) :: eos_yemin  , eos_yemax
     real(kind=c_double) :: energy_shift
     real(kind=c_double) :: temp0, temp1
     real(kind=c_double) :: dlintemp, dlintempi
     real(kind=c_double) :: drholintempi
     real(kind=c_double) :: dlintempyei
     real(kind=c_double) :: drholintempyei
     real(kind=c_double) :: dtemp, dtempi
     real(kind=c_double) :: drho, drhoi
     real(kind=c_double) :: dye, dyei
     real(kind=c_double) :: drhotempi
     real(kind=c_double) :: drhoyei
     real(kind=c_double) :: dtempyei
     real(kind=c_double) :: drhotempyei
  end type NRPyEOS_params
end module NRPyEOS
