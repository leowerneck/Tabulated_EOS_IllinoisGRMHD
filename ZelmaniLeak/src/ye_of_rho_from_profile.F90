module yeofrhomod
  
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  integer nzones

  real*8 rho_min
  real*8 rho_max
  
  real*8,allocatable,save :: yetable(:)
  real*8,allocatable,save :: rhotable(:)
  
  real*8,allocatable,save :: pyea(:)
  real*8,allocatable,save :: prhoa(:)
       
  real*8 density_factor

end module yeofrhomod

subroutine yeofrho(rho,ye)

  use yeofrhomod
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  
  real*8 rho
  real*8 ye
  
  real*8 lrho,buffer
  integer iplus,iminus
  
  integer i

                                                                                                     
! find index in rho/ye tables                                                                        
      lrho = log10(rho)                                                                              
      if(lrho.ge.rhotable(1)) then                                                                   
         ye = yetable(1)                                                                             
         return                                                                                      
      else if (lrho.lt.rhotable(nzones)) then                                                        
         ye = yetable(nzones)                                                                        
         return                                                                                      
      end if                                                                                         
                                                                                                     
      buffer = (lrho-rhotable(1))/(rhotable(nzones)-rhotable(1)) &                                   
           * (nzones-1)*1.0d0                                                                        
      iminus = int(buffer)+1                                                                         
      iplus  = int(buffer)+2                                                                         
                                                                                                     
! call linear interpolation routine                                                                  
      call linterp_ye(rhotable(iplus),rhotable(iminus), &                                            
           yetable(iplus), yetable(iminus),  &                                                       
           lrho, ye )                                                                                

end subroutine yeofrho



subroutine readprofile_ye(profilename,profilezones,prho,pye)

#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  character(*) profilename
  integer profilezones
  integer i

  real*8 bradius
  real*8 prho(*)
  real*8 pye(*)
  
  open(666,file=profilename,status='old')

  do i=1,profilezones
     read(666,*) prho(i),pye(i),bradius
  enddo
  close(666)
  
end subroutine  readprofile_ye

subroutine setuprho_ye(nzones,rho_min,rho_max,rhotable)

#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
! sets up a logarithmically spaced density array                                                               

  real*8 rho_min,rho_max
  real*8 rhotable(*)
  integer nzones

  real*8 lrmin,lrmax,dlr
  integer i

  lrmin = log10(rho_min)
  lrmax = log10(rho_max)

  dlr = (lrmax-lrmin)/((nzones-1)*1.0d0)
  
  rhotable(1) = lrmax
  do i=2,nzones
     rhotable(i) = rhotable(i-1)-dlr
  enddo
  

end subroutine setuprho_ye


subroutine setupye_ye(nzones,profilezones,rhotable,yetable,prho,pye)

#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none
  integer nzones,profilezones
  real*8 rhotable(*),yetable(*),prho(*),pye(*)
  
  real*8 localrho
  
  integer i,upper_index,lower_index

  do i=1,nzones
     localrho=10.0d0**rhotable(i)
     
     if ( localrho.ge.prho(1) ) then
        yetable(i) = pye(1)
     else if(localrho.le.prho(profilezones)) then
        yetable(i) = pye(profilezones)
     else
        !            write(*,"(1P10E15.6)") localrho                                                                   
        call map_find_index_ye(profilezones,prho,localrho, &
             upper_index,lower_index)
        
        call linterp_ye(prho(upper_index),prho(lower_index), &
             pye(upper_index), pye(lower_index),  &
             localrho, yetable(i) )
        
     endif
  enddo

end subroutine setupye_ye


! ******************************************************************                                           
subroutine map_find_index_ye(zones,array,goal,  &
     upper_index,lower_index)

! bisection search                                                                                             
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none

  integer zones,i
  real*8 array(*)
  real*8 goal
  integer middle_index,upper_index,lower_index

  lower_index = 1
  upper_index = zones

  do while ( (upper_index - lower_index) .gt. 1 )
     middle_index = (lower_index + upper_index) * 0.5
     if ( (goal .le. array(lower_index)) &
          .and. (goal .ge. array(middle_index)) ) then
        upper_index = middle_index
     else
        if ( (goal .le. array(middle_index)) &
             .and. (goal .ge. array(upper_index)) ) then
           lower_index = middle_index
        endif
     endif
  enddo
  
end subroutine map_find_index_ye

! ******************************************************************     

subroutine linterp_ye(x1,x2,y1,y2,x,y)
! perform linear interpolation                                                                                 
#ifdef HAVE_CAPABILITY_Fortran
  use cctk
#endif
  implicit none

  real*8 slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     stop "Error in linterp!"
  endif

  slope = (y2 - y1) / (x2 - x1)
  
  !  write(*,*) "slope", slope, x1,x2,y1,y2,x                                                              

  y = slope*(x-x1) + y1

  
end subroutine linterp_ye
