!****************************************************************************
! 磁界に対する吸収境界条件
!****************************************************************************
subroutine absorbing_boundary_condition_for_h
   implicit none

   call mur_1st_for_h        ! Mur 1次の吸収境界条件
   !call mur_2nd_for_h        ! Mur 2次の吸収境界条件

   return
end subroutine

!----------------------------------------------------------------------------
! Mur 1次の吸収境界条件
!----------------------------------------------------------------------------
subroutine mur_1st_for_h
   use consts
   use fdtd
   implicit none

   integer :: i,j,k

   ! **************** 壁 i=1,nx に対して (x軸に沿って) ****************
   ! ---------------- Hy に対して ----------------
   do k=1,nz-1
      do j=2,ny-1
         hy(1,j,k) =hyx1(2,j,k)+cxd*(hy(2,j,k)-hyx1(1,j,k))
         hy(nx-1,j,k)=hyx1(3,j,k)+cxu*(hy(nx-2,j,k)-hyx1(4,j,k))
      end do
   end do
   ! 過去の時間の値の更新
   do k=1,nz-1
      do j=2,ny-1
         hyx1(1,j,k)=hy(1,j,k)
         hyx1(2,j,k)=hy(2,j,k)
         hyx1(3,j,k)=hy(nx-2,j,k)
         hyx1(4,j,k)=hy(nx-1,j,k)
      end do
   end do
   ! ---------------- Hz に対して ----------------
   do k=2,nz-1
      do j=1,ny-1
         hz(1,j,k) =hzx1(2,j,k)+cxd*(hz(2,j,k)-hzx1(1,j,k))
         hz(nx-1,j,k)=hzx1(3,j,k)+cxu*(hz(nx-2,j,k)-hzx1(4,j,k))
      end do
   end do
   ! 過去の時間の値の更新
   do k=2,nz-1
      do j=1,ny-1
         hzx1(1,j,k)=hz(1,j,k)
         hzx1(2,j,k)=hz(2,j,k)
         hzx1(3,j,k)=hz(nx-2,j,k)
         hzx1(4,j,k)=hz(nx-1,j,k)
      end do
   end do

   ! **************** 壁 j=1,ny に対して (y軸に沿って) ****************
   ! ---------------- Hx に対して ----------------
   do k=1,nz-1
      do i=2,nx-1
         hx(i,1,k) =hxy1(i,2,k)+cyd*(hx(i,2,k)-hxy1(i,1,k))
         hx(i,ny-1,k)=hxy1(i,3,k)+cyu*(hx(i,ny-2,k)-hxy1(i,4,k))
      end do
   end do
   ! 過去の時間の値の更新
   do k=1,nz-1
      do i=2,nx-1
         hxy1(i,1,k)=hx(i,1,k)
         hxy1(i,2,k)=hx(i,2,k)
         hxy1(i,3,k)=hx(i,ny-2,k)
         hxy1(i,4,k)=hx(i,ny-1,k)
      end do
   end do
   ! ---------------- Hz に対して ----------------
   do k=2,nz-1
      do i=1,nx-1
         hz(i,1,k) =hzy1(i,2,k)+cyd*(hz(i,2,k)-hzy1(i,1,k))
         hz(i,ny-1,k)=hzy1(i,3,k)+cyu*(hz(i,ny-2,k)-hzy1(i,4,k))
      end do
   end do
   ! 過去の時間の値の更新
   do k=2,nz-1
      do i=1,nx-1
         hzy1(i,1,k)=hz(i,1,k)
         hzy1(i,2,k)=hz(i,2,k)
         hzy1(i,3,k)=hz(i,ny-2,k)
         hzy1(i,4,k)=hz(i,ny-1,k)
      end do
   end do

   ! **************** 壁 k=1,nz に対して (z軸に沿って) ****************
   ! ---------------- Hx に対して ----------------
   do i=2,nx-1
      do j=1,ny-1
         hx(i,j,1) =hxz1(i,j,2)+czd*(hx(i,j,2)-hxz1(i,j,1))
         hx(i,j,nz-1)=hxz1(i,j,3)+czu*(hx(i,j,nz-2)-hxz1(i,j,4))
      end do
   end do
   ! 過去の時間の値の更新
   do i=2,nx-1
      do j=1,ny-1
         hxz1(i,j,1)=hx(i,j,1)
         hxz1(i,j,2)=hx(i,j,2)
         hxz1(i,j,3)=hx(i,j,nz-2)
         hxz1(i,j,4)=hx(i,j,nz-1)
      end do
   end do
   ! ---------------- Hy に対して ----------------
   do i=1,nx-1
      do j=2,ny-1
         hy(i,j,1) =hyz1(i,j,2)+czd*(hy(i,j,2)-hyz1(i,j,1))
         hy(i,j,nz-1)=hyz1(i,j,3)+czu*(hy(i,j,nz-2)-hyz1(i,j,4))
      end do
   end do
   ! 過去の時間の値の更新
   do i=1,nx-1
      do j=2,ny-1
         hyz1(i,j,1)=hy(i,j,1)
         hyz1(i,j,2)=hy(i,j,2)
         hyz1(i,j,3)=hy(i,j,nz-2)
         hyz1(i,j,4)=hy(i,j,nz-1)
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! Mur 2次の吸収境界条件 (p.62)
!----------------------------------------------------------------------------
subroutine mur_2nd_for_h
   implicit none

   !
   ! Mur の2次吸収境界条件 (p.65)
   !
   !      ----------------
   !  z  /|     ④      /|
   !   /  |    (6)    /  |
   !  ----------------   |
   !  |(2)|y         |(5)|
   !  |   -----------|----
   !  |  /    (1)    |  /
   !  |/   ③        |/
   !  ---------------- x
   !
   ! ○: 表の面, (*): 隠れた面
   ! 1: z=0, x-y plane
   ! 2: x=0, y-z plane
   ! 3: y=0, x-z plane
   ! 4: z=z, x-y plane
   ! 5: x=x, y-z plane
   ! 6: y=y, x-z plane
   !

   !call mur_2nd_yz_plane_for_h    ! y-z plane (2,5) 平面に対して
   !call mur_2nd_xz_plane_for_h    ! x-z plane (3,6) 平面に対して
   !call mur_2nd_xy_plane_for_h    ! x-y plane (1,4) 平面に対して

   return
end subroutine




