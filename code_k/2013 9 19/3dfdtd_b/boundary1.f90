!****************************************************************************
! 補正境界条件
!****************************************************************************

!----------------------------------------------------------------------------
! 電界の境界条件
!----------------------------------------------------------------------------
subroutine electric_boundary_condition
   use fdtd
   implicit none

   integer :: i,j,k

   ! 下の壁
   k=1
   do i=1,nx-1
      do j=1,ny
         ex(i,j,k)=0.0d0
      end do
   end do
   do i=1,nx
      do j=1,ny-1
         ey(i,j,k)=0.0d0
      end do
   end do

   ! 上の壁
   k=nz-1
   do i=1,nx-1
      do j=1,ny
         ex(i,j,k)=0.0d0
      end do
   end do
   do i=1,nx
      do j=1,ny-1
         ey(i,j,k)=0.0d0
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! 磁界の境界条件
!----------------------------------------------------------------------------
subroutine magnetic_boundary_condition
   use fdtd
   implicit none

   integer :: i,j,k

   ! 右の壁
   i=nx-1
   do j=1,ny
      do k=1,nz-1
         hy(i,j,k)=0.0d0
      end do
   end do
   do j=1,ny-1
      do k=1,nz
         hz(i,j,k)=0.0d0
      end do
   end do

   ! 左の壁
   i=1
   do j=1,ny
      do k=1,nz-1
         hy(i,j,k)=0.0d0
      end do
   end do
   do j=1,ny-1
      do k=1,nz
         hz(i,j,k)=0.0d0
      end do
   end do

   return
end subroutine

!
! End of file
