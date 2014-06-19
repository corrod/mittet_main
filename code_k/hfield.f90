!-------------------------------------------------------------------------------
! 磁界の計算
!-------------------------------------------------------------------------------
subroutine magnetic_field
  use fdtd
  implicit none

  integer :: i,j,k,id

  ! Hx
  do k=1,nz-1
     do j=1,ny-1
        do i=2,nx-1
           id=media_id(i,j,k)

           if(id.eq.1) then
              ! 1:自由空間
              hx(i,j,k)=hx(i,j,k)-chxry0*(ez(i,j+1,k)-ez(i,j,k))+chxrz0*(ey(i,j,k+1)-ey(i,j,k))
           else if(id.eq.2) then
              ! 2: 完全導体
              hx(i,j,k)=0.0d0
           else
              ! 3以上: 任意媒質
              hx(i,j,k)=hx(i,j,k)-chxry(id)*(ez(i,j+1,k)-ez(i,j,k))+chxrz(id)*(ey(i,j,k+1)-ey(i,j,k))
           end if
        end do
     end do
  end do

  ! Hy
  do k=1,nz-1
     do j=2,ny-1
        do i=1,nx-1
           id=media_id(i,j,k)
           
           if(id.eq.1) then
              ! 1: 自由空間
              hy(i,j,k)=hy(i,j,k)-chyrz0*(ex(i,j,k+1)-ex(i,j,k))+chyrx0*(ez(i+1,j,k)-ez(i,j,k))
           else if(id.eq.2) then
              ! 2: 完全導体
              hy(i,j,k)=0.0d0
           else
              ! 3以上: 任意媒質
              hy(i,j,k)=hy(i,j,k)-chyrz(id)*(ex(i,j,k+1)-ex(i,j,k))+chyrx(id)*(ez(i+1,j,k)-ez(i,j,k))
           end if
        end do
     end do
  end do

! Hz
  do k=2,nz-1
     do j=1,ny-1
        do i=1,nx-1
           id=media_id(i,j,k)

           if(id.eq.1) then
              ! 1: 自由空間
              hz(i,j,k)=hz(i,j,k)-chzry0*(ey(i+1,j,k)-ey(i,j,k))+chzry0*(ex(i,j+1,k)-ex(i,j,k))
           else if(id.eq.2) then
              ! 2: 完全導体
              hz(i,j,k)=0.0d0
           else
              ! 3以上: 任意媒質
              hz(i,j,k)=hz(i,j,k)-chzrx(id)*(ey(i+1,j,k)-ey(i,j,k))+chzry(id)*(ex(i,j+1,k)-ex(i,j,k))
           end if
        end do
     end do
  end do

  return
end subroutine
