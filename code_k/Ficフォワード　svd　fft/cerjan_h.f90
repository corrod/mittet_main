
subroutine cerjan_h 
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


  ! i=1

  do i=2,l
     do j=1,ny-1
        do k=1,nz-1
           hx(i,j,k)=hx(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=1,l
     do j=2,ny-1
        do k=1,nz-1
           hy(i,j,k)=hy(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=1,l
     do j=1,ny-1
        do k=2,nz-1
           hz(i,j,k)=hz(i,j,k)*cerx(i)
        end do
     end do
  end do

  ! i=nx

  do i=nx-l+2,nx-1
     do j=2,ny-1
        do k=2,nz-1
           hx(i,j,k)=hx(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=nx-l+1,nx-1
     do j=1,ny-1
        do k=2,nz-1
           hy(i,j,k)=hy(i,j,k)*cerx(i)
        end do
     end do
  end do
 do i=nx-l+1,nx-1
     do j=2,ny-1
        do k=1,nz-1
           hz(i,j,k)=hz(i,j,k)*cerx(i)
        end do
     end do
  end do


  ! j=1 

  do i=2,nx-1
     do j=1,l
        do k=1,nz-1
           hx(i,j,k)=hx(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=1,nx-1
     do j=2,l
        do k=1,nz-1
           hy(i,j,k)=hy(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=1,nx-1
     do j=1,l
        do k=2,nz-1
           hz(i,j,k)=hz(i,j,k)*cery(j)
        end do
     end do
  end do

  ! j=ny

  do i=2,nx-1
     do j=ny-l+1,ny-1
        do k=1,nz-1
           hx(i,j,k)=hx(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=1,nx-1
     do j=ny-l+2,ny-1
        do k=1,nz-1
           hy(i,j,k)=hy(i,j,k)*cery(j)
        end do
     end do
  end do
 do i=1,nx-1
     do j=ny-l+1,ny-1
        do k=2,nz-1
           hz(i,j,k)=hz(i,j,k)*cery(j)
        end do
     end do
  end do


  ! k=1

  do i=2,nx-1
     do j=1,ny-1
        do k=1,l
           hx(i,j,k)=hx(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=1,nx-1
     do j=2,ny-1
        do k=1,l
           hy(i,j,k)=hy(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=1,nx-1
     do j=1,ny-1
        do k=2,l
           hz(i,j,k)=hz(i,j,k)*cerz(k)
        end do
     end do
  end do

  ! k=nz

  do i=2,nx-1
     do j=1,ny-1
        do k=nz-l+1,nz-1
           hx(i,j,k)=hx(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=1,nx-1
     do j=2,ny-1
        do k=nz-l+1,nz-1
           hy(i,j,k)=hy(i,j,k)*cerz(k)
        end do
     end do
  end do
 do i=1,nx-1
     do j=1,ny-1
        do k=nz-l+2,nz-1
           hz(i,j,k)=hz(i,j,k)*cerz(k)
        end do
     end do
  end do

end subroutine
