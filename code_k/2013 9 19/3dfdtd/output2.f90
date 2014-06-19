!******************************************************************************
! 出力 タイムラプス
!******************************************************************************
subroutine output2
  use fdtd
  implicit none

  integer :: i,j,k

   open (40,file="receiver1.dat")
   i=nx/2
   j=ny/2
   k=nz/4
!  write(40,*) n*dt,sqrt(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
   write(40,*) n,ey(i,j,k)

   return
end subroutine
