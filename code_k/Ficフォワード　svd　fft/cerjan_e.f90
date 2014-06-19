
subroutine cerjan_e 
  use fdtd 
  implicit none 
  integer :: i,j,k,l
  real(8) :: a,cerx(mx),cery(my),cerz(mz)

  do i=1,nx
     cerx(i)=1.0d0
  end do
  do j=1,ny
     cery(j)=1.0d0
  end do
  do k=1,nz
     cerz(k)=1.0d0
  end do

  l=30
  a=5.0d-5
  do i=1,l
     cerx(i)=dexp(-a*(l-i)**2)
     cerx(nx-i+1)=cerx(i)
  end do

  do j=1,l
     cery(j)=dexp(-a*(l-j)**2)
     cery(ny-j+1)=cery(j)
  end do

  do k=1,l
     cerz(k)=dexp(-a*(l-k)**2)
     cerz(nz-k+1)=cerz(k)
  end do

 
  open (11,file='cerx(i)') 
  open (12,file='cery(j)') 
  open (13,file='cerz(k)') 

  do i=1,nx
     write(11,*) i,cerx(i)
  end do
  do j=1,ny
     write(12,*) j,cery(j)
  end do
  do k=1,nz
     write(13,*) k,cerz(k)
  end do
  close(11)
  close(12)
  close(13)

  ! i=1

  do i=1,l
     do j=2,ny-1
        do k=2,nz-1
           ex(i,j,k)=ex(i,j,k)*cerx(i)
           
        end do
     end do
  end do
 do i=2,l
     do j=1,ny-1
        do k=2,nz-1
           ey(i,j,k)=ey(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=2,l
     do j=2,ny-1
        do k=1,nz-1
           ez(i,j,k)=ez(i,j,k)*cerx(i)
        end do
     end do
  end do

  ! i=nx
  do i=nx-l+1,nx-1
     do j=2,ny-1
        do k=2,nz-1
           ex(i,j,k)=ex(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=nx-l+2,nx-1
     do j=1,ny-1
        do k=2,nz-1
           ey(i,j,k)=ey(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=nx-l+2,nx-1
     do j=2,ny-1
        do k=1,nz-1
           ez(i,j,k)=ez(i,j,k)*cerx(i)
        end do
     end do
  end do

  ! j=1 

  do i=1,nx-1
     do j=2,l
        do k=2,nz-1
           ex(i,j,k)=ex(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=2,nx-1
     do j=1,l
        do k=2,nz-1
           ey(i,j,k)=ey(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=2,nx-1
     do j=2,l
        do k=1,nz-1
           ez(i,j,k)=ez(i,j,k)*cery(j)
        end do
     end do
  end do

  ! j=ny
  do i=1,nx-1
     do j=ny-l+2,ny-1
        do k=2,nz-1
           ex(i,j,k)=ex(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=2,nx-1
     do j=ny-l+1,ny-1
        do k=2,nz-1
           ey(i,j,k)=ey(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=2,nx-1
     do j=ny-l+2,ny-1
        do k=1,nz-1
           ez(i,j,k)=ez(i,j,k)*cery(j)
        end do
     end do
  end do

  ! k=1

  do i=1,nx-1
     do j=2,ny-1
        do k=2,l
           ex(i,j,k)=ex(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=2,nx-1
     do j=1,ny-1
        do k=2,l
           ey(i,j,k)=ey(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=2,nx-1
     do j=2,ny-1
        do k=1,l
           ez(i,j,k)=ez(i,j,k)*cerz(k)
        end do
     end do
  end do

  ! k=nz

  do i=1,nx-1
     do j=2,ny-1
        do k=nz-l+2,nz-1
           ex(i,j,k)=ex(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=2,nx-1
     do j=1,ny-1
        do k=nz-l+2,nz-1
           ey(i,j,k)=ey(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=2,nx-1
     do j=2,ny-1
        do k=nz-l+1,nz-1
           ez(i,j,k)=ez(i,j,k)*cerz(k)
        end do
     end do
  end do

end subroutine
