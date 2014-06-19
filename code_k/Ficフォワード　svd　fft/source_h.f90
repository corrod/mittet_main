
!-------------------------------------------------------------------------------
!    
!-------------------------------------------------------------------------------

subroutine ecur_infinitesimal_dipole_source
   use consts 
   use fdtd
   implicit none

   integer :: i,j,k,id
   real(8) :: sin_warming_up, func_gaussian, func_sin, func_ricker,func_gaussian1

   i=nx/2
!   do j=ny/2-18-2,ny/2-18+2
   do j=ny/2-35-2,ny/2-35+2
!   j=ny/2-35
       k=nz/2+16
!      k=nz/2+1
      id=media_id(i,j,k)

      ey(i,j,k)=ey(i,j,k)-(dt/eps_p(id))*(1.0d0/(dx*dy))*func_gaussian1(time-dt/2.0d0)
!      ey(i,j,k)=ey(i,j,k)-(dt/eps_p(id))*sin_warming_up(time-dt/2.0d0)
   end do
   return
end subroutine



!-------------------------------------------------------------------------------
!       
!          
!-------------------------------------------------------------------------------
real(8) function func_gaussian1(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp
  real(8) :: beta

  tp=pi/(freq*12.5d0)
  beta=pi*((freq*12.5d0)**2)

  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)/237.26829917425468

  end function
