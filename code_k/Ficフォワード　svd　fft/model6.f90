!------------------------------------------------------------------------------
!          
!------------------------------------------------------------------------------
subroutine modeling
  use fdtd
  implicit none

  integer :: i,j,k
 
  !             
  do i=1,nx
     do j=1,ny
        do k=1,nz
           media_id(i,j,k)=3
        end do
     end do
  end do

  !    1
    do i=1,nx
     do j=1,ny
        do k=1,nz/2+15
           media_id(i,j,k)=5
        end do
     end do
  end do  

return
end subroutine

!-------------------------------------------------------------------------------!         
!-------------------------------------------------------------------------------
subroutine test_modeling
  use fdtd
  implicit none
 
  integer :: i,j,k

  open(15,file="test_modeling.dat")
  j=ny/2
  do i=1,nx
     do k=1,nz
        write(15,*) i,k,media_id(i,j,k)
     end do
     write(15,*)
  end do
 close(15)

  return
end subroutine
