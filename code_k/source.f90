!****************************************************************************
! 電流源
!****************************************************************************
!----------------------------------------------------------------------------
! 線電荷
!----------------------------------------------------------------------------

subroutine ecur_line_source
   use fdtd
   implicit none

   integer :: i,j,k,id
   real(8) :: sin_warming_up, func_gaussian

   i=nx/2
   j=ny/2
   do k=2,nz-1
      id=media_id(i,j,k)

      ey(i,j,k)=ey(i,j,k) &
               -(dt/eps(id))/(1.0d0+(sig(id)*dt)/(2.0d0*eps(id)))*func_gaussian(time-dt/2.0d0)
   end do

   return
end subroutine

!-------------------------------------------------------------------------------
! 点電荷
!-------------------------------------------------------------------------------

subroutine ecur_infinitesimal_dipole_source
   use consts
   use fdtd
   implicit none

   integer :: i,j,k,id
   real(8) :: sin_warming_up, func_gaussian, func_sin, func_ricker,func_gaussian1

   i=nx/2
   j=ny/2
   k=nz/2

   id=media_id(i,j,k)
   ey(i,j,k)=ey(i,j,k) &
            -(dt/eps(id))/(1.0d0+(sig(id)*dt)/(2.0d0*eps(id)))*func_gaussian1(time-dt/2.0d0)

   return
end subroutine

!----------------------------------------------------------------------------
! 平面波
!----------------------------------------------------------------------------
subroutine plane_wave_source
   use consts
   use fdtd
   implicit none

   integer :: i,j,k
   real(8) :: sin_warming_up, func_sin

   j=3
   do k=1,nz-1
      do i=1,nx-1
         ez(i,j,k)=func_sin(time-dt/2.0d0)
         hx(i,j,k)=ez(i,j,k)/z0
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! 波源の関数形
! 正弦波
!----------------------------------------------------------------------------
real(8) function func_sin(t)
   use consts
   use fdtd
   implicit none

   integer :: j
   real(8) :: t
   real(8) :: tp,t1

   tp=1.0d0/freq
!   t1=4.0d0*tp
!   if(t < t1) then
!      ! ウォーミングアップ
!      sin_warming_up=0.5d0*(1.0d0-dcos(pi*t/t1))*dsin(2.0d0*pi*freq*t)
!   else
      ! 正弦波
      func_sin=-dsin(2.0d0*pi*freq*t)
!   end if

end function

!-------------------------------------------------------------------------------
! 波源の関数形
! だんだん振幅が大きくなる正弦波
!-------------------------------------------------------------------------------
real(8) function sin_warming_up(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp,t1

  tp=1.0d0/freq
  t1=3.0d0*tp
  if(t.lt.t1) then
     sin_warming_up=0.5d0*(1.0d0-dcos(pi*t/t1))*dsin(2.0d0*pi*freq*t)
  else
     sin_warming_up=dsin(2.0d0*pi*freq*t)
  end if
end function

!-------------------------------------------------------------------------------
! 波源の関数形
! 余弦波
!-------------------------------------------------------------------------------
real(8) function func_cos(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t

  func_cos=dcos(2.0d0*pi*freq*t)
end function

!-------------------------------------------------------------------------------
! 波源の関数形
! ステップパルス
!-------------------------------------------------------------------------------
real(8) function func_step(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp

  tp=1.0d0/freq
  if(t.lt.(tp/2.0d0)) then
     func_step=1.0d0
  else
     func_step=0.0d0
  end if
end function

!-------------------------------------------------------------------------------
! 波源の関数形
! ガウシアンパルス
!-------------------------------------------------------------------------------
real(8) function func_gaussian(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp,beta

  tp=pi/freq
  beta=pi*freq**2

     func_gaussian=dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)

end function

!-------------------------------------------------------------------------------
! 波源の関数形
! ガウシアン一次微分
!-------------------------------------------------------------------------------
real(8) function func_gaussian1(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp
  real(8) :: beta

  tp=pi/freq
  beta=pi*freq**2

  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)
  end function

!------------------------------------------------------------------------------
! 波源の関数形
! Ricker wavelet
!------------------------------------------------------------------------------
real(8) function func_ricker(t)
  use consts
  use fdtd
  implicit none

  real(8) :: t
  real(8) :: tp

  tp=1.0d0/freq
!  if (t.lt.tp) then
     func_ricker=(1.0d0-2*(pi*freq*(t-tp/2.0d0))**2)*dexp(-(pi*freq*(t-tp/2.0d0))**2)
!  else
!     func_ricker=0.0d0
!  end if
end function


