#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

integer function EOS_Omni_Startup()

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_PARAMETERS

  integer :: p

  if(poly_gamma_initial .gt. 0d0) then
    poly_gamma_ini = poly_gamma_initial
  else
    poly_gamma_ini = poly_gamma
  end if

  poly_k_cgs = poly_k * rho_gf**poly_gamma_ini / press_gf

  gl_k_cgs   = gl_k * rho_gf**poly_gamma_ini / press_gf

! --------------------------------------------
! Use Legacy parameter in case n_pieces=2 
!       and hybrid_gamma(1) has not been set
! --------------------------------------------
    
  if(n_pieces .gt. 0) then
     allocate(hybrideos_eps(n_pieces))
     allocate(hybrideos_k(n_pieces))
     allocate(hybrideos_gamma(n_pieces))
     allocate(hybrideos_rho(n_pieces))
     if( (n_pieces .eq. 2) .and. (hybrid_gamma(1) .eq. 0d0 ) )then
        ! backwards comaptible behaviour for 2 piece hybrid EOS
        hybrideos_k(1)     = hybrid_k1
        if (poly_gamma_initial .gt. 0d0) then
          hybrideos_k(1)   = hybrideos_k(1) * rho_gf**(hybrideos_gamma(1) - poly_gamma_initial)
        end if
        hybrideos_rho(1)   = hybrid_rho_nuc
        hybrideos_gamma(1) = hybrid_gamma1
        hybrideos_gamma(2) = hybrid_gamma2
     else
        ! regular behaviour for n polytropic pieces
        hybrideos_k(1) = hybrid_k0
        do p = 1,n_pieces
           hybrideos_gamma(p) = hybrid_gamma(p)
           hybrideos_rho(p) = hybrid_rho(p)
        end do
     end if
     hybrideos_eps(1) = 0
     if (n_pieces .gt. 1) then
        do p = 1,n_pieces-1
           hybrideos_k(p+1) = hybrideos_k(p) * hybrideos_rho(p)**(hybrideos_gamma(p)-hybrideos_gamma(p+1))
           hybrideos_eps(p+1) = hybrideos_eps(p) + &
                hybrideos_k(p)*hybrideos_rho(p)**(hybrideos_gamma(p)-1.0d0)/(hybrideos_gamma(p)-1.0d0) - &
                hybrideos_k(p+1) * hybrideos_rho(p)**(hybrideos_gamma(p+1)-1.0d0)/(hybrideos_gamma(p+1)-1.0d0)
        end do
     end if
  end if

  EOS_Omni_Startup = 0

end function EOS_Omni_Startup

subroutine EOS_Omni_Get_Energy_Shift(CCTK_ARGUMENTS)

  use EOS_Omni_Module
  implicit none

  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS_EOS_OMNI_GET_ENERGY_SHIFT

  integer :: nuceos_table_was_read_int

  call nuc_eos_c_get_energy_shift(energy_shift,eos_tempmin,eos_tempmax,&
       eos_yemin,eos_yemax,nuceos_table_was_read_int)

  if(nuceos_table_was_read_int .ne. 0) then
     nuceos_table_was_read = .true.
  else
     nuceos_table_was_read = .false.
  end if

end subroutine EOS_Omni_Get_Energy_Shift
