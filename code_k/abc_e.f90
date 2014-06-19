!-------------------------------------------------------------------------------
! 電界に対する吸収境界条件
!-------------------------------------------------------------------------------
subroutine absorbing_boundary_condition_for_e
  implicit none

!  call mur_1st_for_e    ! Mur 1次の吸収境界条件
  call mur_2nd_for_e      ! Mur 2次の吸収境界条件

  return
end subroutine

!-------------------------------------------------------------------------------
! Mur 1次の吸収境界条件
!-------------------------------------------------------------------------------
subroutine mur_1st_for_e
  use consts
  use fdtd
  implicit none

  integer :: i,j,k

  ! *** 壁 i=1,nx に対して（x軸に沿って） ***
  ! --- Eyに対して ---
  do k=2,nz-1
     do j=1,ny-1
        ey(1,j,k) =eyx1(2,j,k)+cxd*(ey(2,j,k)-eyx1(1,j,k))
        ey(nx,j,k)=eyx1(3,j,k)+cxu*(ey(nx-1,j,k)-eyx1(4,j,k))
     end do
  end do
  ! 過去の時間の値の更新
  do k=2,nz-1
     do j=1,ny-1
        eyx1(1,j,k)=ey(1,j,k)
        eyx1(2,j,k)=ey(2,j,k)
        eyx1(3,j,k)=ey(nx-1,j,k)
        eyx1(4,j,k)=ey(nx,j,k)
     end do
  end do

  ! --- Ezに対して ---

  do k=1,nz-1
     do j=2,ny-1
        ez(1,j,k) =ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
        ez(nx,j,k)=ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
     end do
  end do
  ! 過去の時間の値の更新
  do k=1,nz-1
     do j=2,ny-1
        ezx1(1,j,k)=ez(1,j,k)
        ezx1(2,j,k)=ez(2,j,k)
        ezx1(3,j,k)=ez(nx-1,j,k)
        ezx1(4,j,k)=ez(nx,j,k)
     end do
  end do

  ! *** 壁 j=1,ny に対して（y軸に沿って） ***
  ! --- Ex に対して ---
  do k=2,nz-1
     do i=1,nx-1
        ex(i,1,k)=exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
        ex(i,ny,k)=exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))
     end do
  end do

  ! --- Ez に対して ---
  do k=1,nz-1
     do i=2,nx-1
        ez(i,1,k)=ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
        ez(1,ny,k)=ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))
     end do
  end do
  ! 過去の時間の値の更新
  do k=1,nz-1
     do i=2,nx-1
        ezy1(i,1,k)=ez(i,1,k)
        ezy1(i,2,k)=ez(i,2,k)
        ezy1(i,3,k)=ez(i,ny-1,k)
        ezy1(i,4,k)=ez(i,ny,k)
     end do
  end do

  ! *** 壁 k=1,nz に対して（z軸に沿って） ***
  ! --- Ex に対して ---
  do i=1,nx-1
     do j=2,ny-1
        ex(i,j,1)=exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
        ex(i,j,nz)=exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))
     end do
  end do
  ! 過去の時間の値の更新
  do i=1,nx-1
     do j=2,ny-1
        exz1(i,j,1)=ex(i,j,1)
        exz1(i,j,2)=ex(i,j,2)
        exz1(i,j,3)=ex(i,j,nz-1)
        exz1(i,j,4)=ex(i,j,nz)
     end do
  end do
  ! --- Ey に対して ---
  do i=2,nx-1
     do j=1,ny-1
        ey(i,j,1)=eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
        ey(i,j,nz)=eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))
     end do
  end do
  ! 過去の時間の値の更新
  do i=2,nx-1
     do j=1,ny-1
        eyz1(i,j,1)=ey(i,j,1)
        eyz1(i,j,2)=ey(i,j,2)
        eyz1(i,j,3)=ey(i,j,nz-1)
        eyz1(i,j,4)=ey(i,j,nz)
     end do
  end do

return
end subroutine

!-------------------------------------------------------------------------------
!  Mur 2次の吸収境界条件
!-------------------------------------------------------------------------------
subroutine mur_2nd_for_e
  implicit none

  call mur_2nd_yz_plane_for_e         ! y-z plane に対して
  call mur_2nd_xz_plane_for_e         ! x-z plane に対して
  call mur_2nd_xy_plane_for_e         ! x-y plane に対して

  return
end subroutine

!-------------------------------------------------------------------------------
! y-z 平面に対して（x方向に伝播）
!-------------------------------------------------------------------------------
subroutine mur_2nd_yz_plane_for_e
  use consts
  use fdtd
  implicit none

  integer :: i,j,k

  ! *** 壁 i=1,nx に対して ***
  ! --- Ey に対して ---
  ! 1次吸収境界条件
  ! 紫
  do k=2,nz-1
     j=1
     ey(1,j,k) =eyx1(2,j,k)+cyd*(ey(2,j,k)-eyx1(1,j,k))
     ey(nx,j,k)=eyx1(3,j,k)+cyu*(ey(nx-1,j,k)-eyx1(4,j,k))

     j=ny-1
     ey(1,j,k) =eyx1(2,j,k)+cyd*(ey(2,j,k)-eyx1(1,j,k))
     ey(nx,j,k)=eyx1(3,j,k)+cyu*(ey(nx-1,j,k)-eyx1(4,j,k))
  end do
  ! 水色
  do j=2,ny-2
     k=2
     ey(1,j,k) =eyx1(2,j,k)+cyd*(ey(2,j,k)-eyx1(1,j,k))
     ey(nx,j,k)=eyx1(3,j,k)+cyu*(ey(nx-1,j,k)-eyx1(4,j,k))

     k=nz-1
     ey(1,j,k) =eyx1(2,j,k)+cyd*(ey(2,j,k)-eyx1(1,j,k))
     ey(nx,j,k)=eyx1(3,j,k)+cyu*(ey(nx-1,j,k)-eyx1(4,j,k))
  end do
  ! 2次吸収境界条件
  ! 緑
  do k=3,nz-2
     do j=2,ny-2
        ey(1,j,k)=-eyx2(2,j,k)+cyd*(ey(2,j,k)+eyx2(1,j,k)) &
                  +cyy*(eyx1(1,j,k)+eyx1(2,j,k)) &
                  +cxfyd*(eyx1(1,j+1,k)-2.0d0*eyx1(1,j,k)+eyx1(1,j-1,k) &
                         +eyx1(2,j+1,k)-2.0d0*eyx1(2,j,k)+eyx1(2,j-1,k)) &
                  +cxfzd*(eyx1(1,j,k+1)-2.0d0*eyx1(1,j,k)+eyx1(1,j,k-1) &
                         +eyx1(2,j,k+1)-2.0d0*eyx1(2,j,k)+eyx1(2,j,k-1))
        ey(nx,j,k)=-eyx2(3,j,k)+cyd*(ey(nx-1,j,k)+eyx2(4,j,k)) &
                   +cyy*(eyx1(4,j,k)+eyx1(3,j,k))&
                   +cxfyd*(eyx1(4,j+1,k)-2.0d0*eyx1(4,j,k)+eyx1(4,j-1,k) &
                          +eyx1(3,j+1,k)-2.0d0*eyx1(3,j,k)+eyx1(3,j-1,k))&
                   +cxfzd*(eyx1(4,j,k+1)-2.0d0*eyx1(4,j,k)+eyx1(4,j,k-1) &
                          +eyx1(3,j,k+1)-2.0d0*eyx1(3,j,k)+eyx1(3,j,k-1))
     end do
  end do

  ! 過去の時間の値の更新
  do k=2,nz-1
     do j=1,ny-1
        eyx2(1,j,k)=eyx1(1,j,k)
        eyx2(2,j,k)=eyx1(2,j,k)
        eyx2(3,j,k)=eyx1(3,j,k)
        eyx2(4,j,k)=eyx1(4,j,k)

        eyx1(1,j,k)=ey(1,j,k)
        eyx1(2,j,k)=ey(2,j,k)
        eyx1(3,j,k)=ey(nx-1,j,k)
        eyx1(4,j,k)=ey(nx,j,k)
     end do
  end do

  ! --- Ez に対して ---
  ! 1次吸収境界条件
  ! 紫

   do k=1,nz-1
      j=2
      ez(1,j,k) =ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
      ez(nx,j,k)=ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))

      j=ny-1
      ez(1,j,k) =ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
      ez(nx,j,k)=ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
   end do
   ! 水色
   do j=3,ny-2
      k=1
      ez(1,j,k) =ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
      ez(nx,j,k)=ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))

      k=nz-1
      ez(1,j,k) =ezx1(2,j,k)+cxd*(ez(2,j,k)-ezx1(1,j,k))
      ez(nx,j,k)=ezx1(3,j,k)+cxu*(ez(nx-1,j,k)-ezx1(4,j,k))
   end do
   ! ２次吸収境界条件
   ! 緑
   do k=2,nz-2
      do j=3,ny-2
         ez(1,j,k) =-ezx2(2,j,k)+cxd*(ez(2,j,k)+ezx2(1,j,k)) &
                    +cxx*(ezx1(1,j,k)+ezx1(2,j,k)) &
                    +cxfyd*(ezx1(1,j+1,k)-2.0d0*ezx1(1,j,k)+ezx1(1,j-1,k) &
                           +ezx1(2,j+1,k)-2.0d0*ezx1(2,j,k)+ezx1(2,j-1,k)) &
                    +cxfzd*(ezx1(1,j,k+1)-2.0d0*ezx1(1,j,k)+ezx1(1,j,k-1) &
                           +ezx1(2,j,k+1)-2.0d0*ezx1(2,j,k)+ezx1(2,j,k-1))
         ez(nx,j,k)=-ezx2(3,j,k)+cxd*(ez(nx-1,j,k)+ezx2(4,j,k)) &
                    +cxx*(ezx1(4,j,k)+ezx1(3,j,k)) &
                    +cxfyd*(ezx1(4,j+1,k)-2.0d0*ezx1(4,j,k)+ezx1(4,j-1,k) &
                           +ezx1(3,j+1,k)-2.0d0*ezx1(3,j,k)+ezx1(3,j-1,k)) &
                    +cxfzd*(ezx1(4,j,k+1)-2.0d0*ezx1(4,j,k)+ezx1(4,j,k-1) &
                           +ezx1(3,j,k+1)-2.0d0*ezx1(3,j,k)+ezx1(3,j,k-1))
      end do
   end do

   ! 過去の時間の値の更新
   do k=1,nz-1
      do j=2,ny-1
         ezx2(1,j,k)=ezx1(1,j,k)
         ezx2(2,j,k)=ezx1(2,j,k)
         ezx2(3,j,k)=ezx1(3,j,k)
         ezx2(4,j,k)=ezx1(4,j,k)

         ezx1(1,j,k)=ez(1,j,k)
         ezx1(2,j,k)=ez(2,j,k)
         ezx1(3,j,k)=ez(nx-1,j,k)
         ezx1(4,j,k)=ez(nx,j,k)
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! x-z 平面に対して(y軸方向に伝播)
!----------------------------------------------------------------------------
subroutine mur_2nd_xz_plane_for_e
   use consts
   use fdtd
   implicit none

   integer :: i,j,k

   ! **************** 壁 j=1,ny に対して ****************
   ! ---------------- Ex に対して ----------------
   ! １次吸収境界条件
   ! 紫
   do k=2,nz-1
      i=1
      ex(i,1,k) =exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
      ex(i,ny,k)=exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))

      i=nx-1
      ex(i,1,k) =exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
      ex(i,ny,k)=exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))
   end do
   ! 水色
   do i=2,nx-2
      k=2
      ex(i,1,k) =exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
      ex(i,ny,k)=exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))

      k=nz-1
      ex(i,1,k) =exy1(i,2,k)+cyd*(ex(i,2,k)-exy1(i,1,k))
      ex(i,ny,k)=exy1(i,3,k)+cyu*(ex(i,ny-1,k)-exy1(i,4,k))
   end do
   ! ２次吸収境界条件
   ! 緑
   do k=3,nz-2
      do i=2,nx-2
         ex(i,1,k) =-exy2(i,2,k)+cyd*(ex(i,2,k)+exy2(i,1,k)) &
                    +cyy*(exy1(i,1,k)+exy1(i,2,k)) &
                    +cyfxd*(exy1(i+1,1,k)-2.0d0*exy1(i,1,k)+exy1(i-1,1,k) &
                           +exy1(i+1,2,k)-2.0d0*exy1(i,2,k)+exy1(i-1,2,k)) &
                    +cyfzd*(exy1(i,1,k+1)-2.0d0*exy1(i,1,k)+exy1(i,1,k-1) &
                           +exy1(i,2,k+1)-2.0d0*exy1(i,2,k)+exy1(i,2,k-1))
         ex(i,ny,k)=-exy2(i,3,k)+cyd*(ex(i,ny-1,k)+exy2(i,4,k)) &
                    +cyy*(exy1(i,4,k)+exy1(i,3,k)) &
                    +cyfxd*(exy1(i+1,4,k)-2.0d0*exy1(i,4,k)+exy1(i-1,4,k) &
                           +exy1(i+1,3,k)-2.0d0*exy1(i,3,k)+exy1(i-1,3,k)) &
                    +cyfzd*(exy1(i,4,k+1)-2.0d0*exy1(i,4,k)+exy1(i,4,k-1) &
                           +exy1(i,3,k+1)-2.0d0*exy1(i,3,k)+exy1(i,3,k-1))
      end do
   end do

   ! 過去の時間の値の更新
   do k=2,nz-1
      do i=1,nx-1
         exy2(i,1,k)=exy1(i,1,k)
         exy2(i,2,k)=exy1(i,2,k)
         exy2(i,3,k)=exy1(i,3,k)
         exy2(i,4,k)=exy1(i,4,k)

         exy1(i,1,k)=ex(i,1,k)
         exy1(i,2,k)=ex(i,2,k)
         exy1(i,3,k)=ex(i,ny-1,k)
         exy1(i,4,k)=ex(i,ny,k)
      end do
   end do

   ! ---------------- Ez に対して ----------------
   ! １次吸収境界条件
   ! 紫
   do k=1,nz-1
      i=2
      ez(i,1,k) =ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
      ez(i,ny,k)=ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))

      i=nx-1
      ez(i,1,k) =ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
      ez(i,ny,k)=ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))
   end do
   ! 水色
   do i=3,nx-2
      k=1
      ez(i,1,k) =ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
      ez(i,ny,k)=ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))

      k=nz-1
      ez(i,1,k) =ezy1(i,2,k)+cyd*(ez(i,2,k)-ezy1(i,1,k))
      ez(i,ny,k)=ezy1(i,3,k)+cyu*(ez(i,ny-1,k)-ezy1(i,4,k))
   end do
   ! ２次吸収境界条件
   ! 緑
   do k=2,nz-2
      do i=3,nx-2
         ez(i,1,k) =-ezy2(i,2,k)+cyd*(ez(i,2,k)+ezy2(i,1,k)) &
                    +cyy*(ezy1(i,1,k)+ezy1(i,2,k)) &
                    +cyfxd*(ezy1(i+1,1,k)-2.0d0*ezy1(i,1,k)+ezy1(i-1,1,k) &
                           +ezy1(i+1,2,k)-2.0d0*ezy1(i,2,k)+ezy1(i-1,2,k)) &
                    +cyfzd*(ezy1(i,1,k+1)-2.0d0*ezy1(i,1,k)+ezy1(i,1,k-1) &
                           +ezy1(i,2,k+1)-2.0d0*ezy1(i,2,k)+ezy1(i,2,k-1))
         ez(i,ny,k)=-ezy2(i,3,k)+cyd*(ez(i,ny-1,k)+ezy2(i,4,k)) &
                    +cyy*(ezy1(i,4,k)+ezy1(i,3,k)) &
                    +cyfxd*(ezy1(i+1,4,k)-2.0d0*ezy1(i,4,k)+ezy1(i-1,4,k) &
                           +ezy1(i+1,3,k)-2.0d0*ezy1(i,3,k)+ezy1(i-1,3,k)) &
                    +cyfzd*(ezy1(i,4,k+1)-2.0d0*ezy1(i,4,k)+ezy1(i,4,k-1) &
                           +ezy1(i,3,k+1)-2.0d0*ezy1(i,3,k)+ezy1(i,3,k-1))
      end do
   end do

   ! 過去の時間の値の更新
   do k=1,nz-1
      do i=2,nx-1
         ezy2(i,1,k)=ezy1(i,1,k)
         ezy2(i,2,k)=ezy1(i,2,k)
         ezy2(i,3,k)=ezy1(i,3,k)
         ezy2(i,4,k)=ezy1(i,4,k)

         ezy1(i,1,k)=ez(i,1,k)
         ezy1(i,2,k)=ez(i,2,k)
         ezy1(i,3,k)=ez(i,ny-1,k)
         ezy1(i,4,k)=ez(i,ny,k)
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! x-y 平面に対して(z軸方向に伝播)
!----------------------------------------------------------------------------
subroutine mur_2nd_xy_plane_for_e
   use consts
   use fdtd
   implicit none

   integer :: i,j,k

   ! **************** 壁 k=1,nz に対して ****************
   ! ---------------- Ex に対して ----------------
   ! １次吸収境界条件
   ! 紫
   do j=2,ny-1
      i=1
      ex(i,j,1) =exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
      ex(i,j,nz)=exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))

      i=nx-1
      ex(i,j,1) =exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
      ex(i,j,nz)=exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))
   end do
   ! 水色
   do i=2,nx-2
      j=2
      ex(i,j,1) =exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
      ex(i,j,nz)=exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))

      j=ny-1
      ex(i,j,1) =exz1(i,j,2)+czd*(ex(i,j,2)-exz1(i,j,1))
      ex(i,j,nz)=exz1(i,j,3)+czu*(ex(i,j,nz-1)-exz1(i,j,4))
   end do
   ! ２次吸収境界条件
   ! 緑
   do i=2,nx-2
      do j=3,ny-2
         ex(i,j,1) =-exz2(i,j,2)+czd*(ex(i,j,2)+exz2(i,j,1)) &
                    +czz*(exz1(i,j,1)+exz1(i,j,2)) &
                    +czfxd*(exz1(i+1,j,1)-2.0d0*exz1(i,j,1)+exz1(i-1,j,1) &
                           +exz1(i+1,j,2)-2.0d0*exz1(i,j,2)+exz1(i-1,j,2)) &
                    +czfyd*(exz1(i,j+1,1)-2.0d0*exz1(i,j,1)+exz1(i,j-1,1) &
                           +exz1(i,j+1,2)-2.0d0*exz1(i,j,2)+exz1(i,j-1,2))
         ex(i,j,nz)=-exz2(i,j,3)+czd*(ex(i,j,nz-1)+exz2(i,j,4)) &
                    +czz*(exz1(i,j,4)+exz1(i,j,3)) &
                    +czfxd*(exz1(i+1,j,4)-2.0d0*exz1(i,j,4)+exz1(i-1,j,4) &
                           +exz1(i+1,j,3)-2.0d0*exz1(i,j,3)+exz1(i-1,j,3)) &
                    +czfyd*(exz1(i,j+1,4)-2.0d0*exz1(i,j,4)+exz1(i,j-1,4) &
                           +exz1(i,j+1,3)-2.0d0*exz1(i,j,3)+exz1(i,j-1,3))
      end do
   end do

   ! 過去の時間の値の更新
   do i=1,nx-1
      do j=2,ny-1
         exz2(i,j,1)=exz1(i,j,1)
         exz2(i,j,2)=exz1(i,j,2)
         exz2(i,j,3)=exz1(i,j,3)
         exz2(i,j,4)=exz1(i,j,4)

         exz1(i,j,1)=ex(i,j,1)
         exz1(i,j,2)=ex(i,j,2)
         exz1(i,j,3)=ex(i,j,nz-1)
         exz1(i,j,4)=ex(i,j,nz)
      end do
   end do

   ! ---------------- Ey に対して ----------------
   ! １次吸収境界条件
   ! 紫
   do j=1,ny-1
      i=2
      ey(i,j,1) =eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
      ey(i,j,nz)=eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))

      i=nx-1
      ey(i,j,1) =eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
      ey(i,j,nz)=eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))
   end do
   ! 水色
   do i=3,nx-2
      j=1
      ey(i,j,1) =eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
      ey(i,j,nz)=eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))

      j=ny-1
      ey(i,j,1) =eyz1(i,j,2)+czd*(ey(i,j,2)-eyz1(i,j,1))
      ey(i,j,nz)=eyz1(i,j,3)+czu*(ey(i,j,nz-1)-eyz1(i,j,4))
   end do
   ! ２次吸収境界条件
   ! 緑
   do i=3,nx-2
      do j=2,ny-2
         ey(i,j,1) =-eyz2(i,j,2)+czd*(ey(i,j,2)+eyz2(i,j,1)) &
                    +czz*(eyz1(i,j,1)+eyz1(i,j,2)) &
                    +czfxd*(eyz1(i+1,j,1)-2.0d0*eyz1(i,j,1)+eyz1(i-1,j,1) &
                           +eyz1(i+1,j,2)-2.0d0*eyz1(i,j,2)+eyz1(i-1,j,2)) &
                    +czfyd*(eyz1(i,j+1,1)-2.0d0*eyz1(i,j,1)+eyz1(i,j-1,1) &
                           +eyz1(i,j+1,2)-2.0d0*eyz1(i,j,2)+eyz1(i,j-1,2))
         ey(i,j,nz)=-eyz2(i,j,3)+czd*(ey(i,j,nz-1)+eyz2(i,j,4)) &
                    +czz*(eyz1(i,j,4)+eyz1(i,j,3)) &
                    +czfxd*(eyz1(i+1,j,4)-2.0d0*eyz1(i,j,4)+eyz1(i-1,j,4) &
                           +eyz1(i+1,j,3)-2.0d0*eyz1(i,j,3)+eyz1(i-1,j,3)) &
                    +czfyd*(eyz1(i,j+1,4)-2.0d0*eyz1(i,j,4)+eyz1(i,j-1,4) &
                           +eyz1(i,j+1,3)-2.0d0*eyz1(i,j,3)+eyz1(i,j-1,3))
      end do
   end do

   ! 過去の時間の値の更新
   do i=2,nx-1
      do j=1,ny-1
         eyz2(i,j,1)=eyz1(i,j,1)
         eyz2(i,j,2)=eyz1(i,j,2)
         eyz2(i,j,3)=eyz1(i,j,3)
         eyz2(i,j,4)=eyz1(i,j,4)

         eyz1(i,j,1)=ey(i,j,1)
         eyz1(i,j,2)=ey(i,j,2)
         eyz1(i,j,3)=ey(i,j,nz-1)
         eyz1(i,j,4)=ey(i,j,nz)
      end do
   end do

   return
end subroutine


