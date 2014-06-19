!-------------------------------------------------------------------------------
! 電界の計算
!-------------------------------------------------------------------------------
subroutine electric_field
  use fdtd
  implicit none

  integer :: i,j,k,id

  ! Ex
  do k=2,nz-1
     do j=2,ny-1
        do i=1,nx-1
           id=media_id(i,j,k)

           if(id.eq.1) then
              ! 1:自由空間
              ex(i,j,k)=cex0*ex(i,j,k)+cexry0*(hz(i,j,k)-hz(i,j-1,k))-cexrz0*(hy(i,j,k)-hy(i,j,k-1))
           else if(id.eq.2) then
              ! 2: 完全導体
              ex(i,j,k)=0.0d0
           else
              ! 3以上: 任意媒質
!              ex(i,j,k)=cex(id)*ex(i,j,k)+cexry(id)*(hz(i,j,k)-hz(i,j-1,k))-cexrz(id)*(hy(i,j,k)-hy(i,j,k-1))
              ex(i,j,k)=ex(i,j,k)+cexry_p(id)*(hz(i,j,k)-hz(i,j-1,k))-cexrz_p(id)*(hy(i,j,k)-hy(i,j,k-1))
           end if
        end do
     end do
  end do

  ! Ey
  do k=2,nz-1
     do j=1,ny-1
        do i=2,nx-1
           id=media_id(i,j,k)

           if(id.eq.1) then
              ! 1:自由空間
              ey(i,j,k)=cey0*ey(i,j,k)+ceyrz0*(hx(i,j,k)-hx(i,j,k-1))-ceyrx0*(hz(i,j,k)-hz(i-1,j,k))
           else if(id.eq.2) then
              ! 2: 完全導体
              ey(i,j,k)=0.0d0
           else
              ! 3以上: 任意媒質
!              ey(i,j,k)=cey(id)*ey(i,j,k)+ceyrz(id)*(hx(i,j,k)-hx(i,j,k-1))-ceyrx(id)*(hz(i,j,k)-hz(i-1,j,k))
              ey(i,j,k)=ey(i,j,k)+ceyrz_p(id)*(hx(i,j,k)-hx(i,j,k-1))-ceyrx_p(id)*(hz(i,j,k)-hz(i-1,j,k))
           end if
        end do
     end do
  end do

  ! Ez
  do k=1,nz-1
     do j=2,ny-1
        do i=2,nx-1
           id=media_id(i,j,k)

           if(id.eq.1) then
           ! 1: 自由空間
           ez(i,j,k)=cez0*ez(i,j,k)+cezrx0*(hy(i,j,k)-hy(i-1,j,k))-cezry0*(hx(i,j,k)-hx(i,j-1,k))
           else if(id.eq.2) then
              ! 2: 完全導体
              ez(i,j,k)=0.0d0
           else 
              ! 3以上: 任意媒質
!              ez(i,j,k)=cez(id)*ez(i,j,k)+cezrx(id)*(hy(i,j,k)-hy(i-1,j,k))-cezry(id)*(hx(i,j,k)-hx(i,j-1,k))
              ez(i,j,k)=ez(i,j,k)+cezrx_p(id)*(hy(i,j,k)-hy(i-1,j,k))-cezry_p(id)*(hx(i,j,k)-hx(i,j-1,k))
           end if
        end do
     end do
  end do

  return
end subroutine
