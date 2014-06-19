!------------------------------------------------------------------------------
! 物体のモデリング
!------------------------------------------------------------------------------
subroutine modeling
  use fdtd
  implicit none

  integer :: i,j,k
 
  ! すべて真空（＝１）に初期化する
  do i=1,nx
     do j=1,ny
        do k=1,nz
           media_id(i,j,k)=3
        end do
     end do
  end do
  
return
end subroutine

!-------------------------------------------------------------------------------
! モデリングの確認
!-------------------------------------------------------------------------------
subroutine test_modeling
  use fdtd
  implicit none
 
  integer :: i,j,k

  open(15,file="test_modeling.dat")
  k=nz/2
  do j=1,ny
     write(10,*) (media_id(i,j,k),i=1,nx)
!     write(10,'()')
!     write(10,'()')
  end do
  close(15)

  return
end subroutine
