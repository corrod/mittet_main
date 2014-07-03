program loop
	implicit none
	integer i,j,k
	real(8) :: a(100,100,100)

	a(1:100,1:100,1:100) =199.0d0
! 	do k=50,60
! 		a(1:100,1:100,k) = 3.0d0
! 	enddo

		a(1:100,1:100,50:60) = 3.0d0

	write(*,*) a(33,33,55), a(22,22,22)
end program loop