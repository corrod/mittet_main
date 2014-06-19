!****************************************************************************
! モデリング
!****************************************************************************

!----------------------------------------------------------------------------
! 物体のモデリング
!----------------------------------------------------------------------------
subroutine modeling
   use fdtd
   implicit none

   integer :: i,j,k

   ! すべて真空 (=1) に初期化する
   do i=1,nx
      do j=1,ny
         do k=1,nz
            media_id(i,j,k)=1
         end do
      end do
   end do

   return
end subroutine

!----------------------------------------------------------------------------
! モデリングの確認
!----------------------------------------------------------------------------
subroutine test_modeling
   use fdtd
   implicit none

   integer :: i,j,k

   open(10,file="test_modeling.dat")
   k=nz/2
   do j=1,ny
      write(10,*) (media_id(i,j,k),i=1,nx)
      write(10,'()')    ! 改行
      write(10,'()')    ! 改行
   end do
   close(10)

   return
end subroutine

!
! End of file
!

