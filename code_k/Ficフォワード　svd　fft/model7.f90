!------------------------------------------------------------------------------
!          
!------------------------------------------------------------------------------
subroutine modeling
  use fdtd
  implicit none

  integer :: i,j,k
             
  do i=1,nx
     do j=1,ny
        do k=1,nz
           media_id(i,j,k)=3
        end do
     end do
  end do

    do i=1,nx
     do j=1,ny
        do k=1,nz/2+15
           media_id(i,j,k)=5
        end do
     end do
  end do

    do i=nx/2-4,nx/2+5
     do j=ny/2-4,ny/2+5
!        do k=3*nz/5-2,3*nz/5
        k=nz/2+11
          media_id(i,j,k)=7
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
