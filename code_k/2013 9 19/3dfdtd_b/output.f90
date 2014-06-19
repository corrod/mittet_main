!****************************************************************************
! 出力
!****************************************************************************
subroutine output
   use fdtd
   implicit none

   integer :: i,j,k

   open(10,file="fieldxy.dat")
   k=nz/2  
      do j=2,ny-1
         do i=2,nx-1
            write(10,*) i,j,ey(i,j,k)
!           write(10,*) i,j,sqrt(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
!      write(10,"()")    ! 改行
         end do
         write(10,"()")    ! 改行
      end do

   open(20,file="fieldyz.dat")
   i=nx/2
   do j=2,ny-1
      do k=2,nz-1
         write(20,*) j,k,ey(i,j,k)
!         write(20,*) j,k,sqrt(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
      end do
      write(20,*)
   end do

   open(30,file="fieldzx.dat")
   j=ny/2
   do i=2,nx-1
      do k=2,nz-1
          write(30,*) i,k,ey(i,j,k)
!        write(30,*) i,k,sqrt(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
      end do
      write(30,*)
   end do

   open (40,file="receiver1.dat")
   i=nx/2
   j=ny/2
   k=nz/4
!  write(40,*) n*dt,sqrt(ex(i,j,k)**2+ey(i,j,k)**2+ez(i,j,k)**2)
   write(40,*) n,ey(i,j,k)

!      write(10,*)
!   end do
!   close(10)

!   open(10,file="field.dat")
!   k=nz/2
!   do j=2,ny-1
!      write(10,*) (dabs(ez(i,j,k))**2,i=2,nx-1)
!      write(10,"()")    ! 改行
!      write(10,"()")    ! 改行
!   end do
!   close(10)

   return
end subroutine
