#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

subroutine EOS_Omni_Check_EOS_Params(eoskey)

  use EOS_Omni_Module
  implicit none
  DECLARE_CCTK_PARAMETERS

  CCTK_INT, INTENT(IN) :: eoskey

  CCTK_REAL, parameter :: zero = 0

  integer :: p
  character(256) :: warnstring

  select case (eoskey)
    case (1)
       ! polytropic EOS
       if(poly_k_cgs .le. zero .or. poly_gamma .le. zero) then
         write(warnstring,*) "invalid values for poly_k_cgs ",&
                             poly_k_cgs, " or poly_gamma ", poly_gamma
         call CCTK_ERROR(warnstring)
         STOP
       endif
    case (2)
       ! gamma-law EOS
       if(gl_k_cgs .le. zero .or. gl_gamma .le. zero) then
         write(warnstring,*) "invalid values for gl_k_cgs ",&
                             gl_k_cgs, " or gl_gamma ", gl_gamma
         call CCTK_ERROR(warnstring)
         STOP
       endif
    case (3)
       ! hybrid EOS
       if(n_pieces .lt. 1) then
         write(warnstring,*) "hybrid EOS needs at least 1 piece, but ",&
                             n_pieces, " were given"
         call CCTK_ERROR(warnstring)
         STOP
       endif
       do p = 1,n_pieces
         if(p .lt. n_pieces) then
           if(hybrideos_rho(p) .le. zero) then
             write(warnstring,*) &
               "separation points must be positive, but point ",p-1,&
               " has value ",hybrideos_rho(p)
             call CCTK_ERROR(warnstring)
           endif
           if(p .gt. 1) then
             if(hybrideos_rho(p) .le. hybrideos_rho(p-1)) then
               write(warnstring,*) &
                 "separation points must be strictly increasing, but hybrid_rho[",&
                 p-1,"] <= hybrid_rho[",p-2,"]: ", hybrideos_rho(p),&
                 " and ",hybrideos_rho(p-1)
               call CCTK_ERROR(warnstring)
               STOP
             endif
           endif
         endif
         ! TODO: is the a constraint on the eps?
         if(hybrideos_gamma(p) .le. zero .or. hybrideos_k(p) .le. zero) then
           write(warnstring, *) "hydrid EOS hybrid_gamma[] and hybrid_kp[] must be strictly positive, but in the set ",&
                                p-1, " they are ", hybrideos_gamma(p)," and ",&
                                hybrideos_k(p)
           call CCTK_ERROR(warnstring)
           STOP
         endif
       end do
    case (4)
       ! nuc eos
       if(.not. nuceos_table_was_read) then
         write(warnstring, *)  "Cold EOS table was not read. Likely you need to set coldeos_read_table"
           call CCTK_ERROR(warnstring)
           STOP
       endif
    case (5)
       ! cold tabular EOS with gamma law
       if(.not. coldeos_table_was_read) then
         write(warnstring, *)  "Cold EOS table was not read. Likely you need to set coldeos_read_table"
           call CCTK_ERROR(warnstring)
           STOP
       endif
    case (6)
       ! barotropic tabular EOS with gamma law
       if(.not. barotropiceos_table_was_read) then
         write(warnstring, *)  "Barotropic EOS table was not read. Likely you need to set barotropiceos_read_table"
           call CCTK_ERROR(warnstring)
           STOP
       endif
    case DEFAULT
       ! nothing to check
       continue
  end select
end subroutine EOS_Omni_Check_EOS_Params
