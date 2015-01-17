program main
	implicit none
	integer i,j

	real(8) ::a(3)
	 a(1)=3.0d0
	 a(2)=4.0d0
	 a(3)=5.0d0

do j=1,100
	open(10,file='outtest.dat',access='append')
	do i = 1,3
		write(10,fmt='(e22.15)', advance='no')  a(i)
    enddo
        close(10)
enddo

write(*,*) 9.0d0*4.5d-3/1767.77d0/1.25846d-6


end program main
! program advance_number
!   implicit none
!   integer,parameter :: n = 7
!   integer i
!   do i = 1, n
!      write (*, fmt='(I4)', advance='no') i   ! 改行しない
!   end do
! end program advance_number