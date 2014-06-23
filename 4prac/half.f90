program half
	implicit none
	integer :: i,j,k
	integer ,parameter :: L=10, M=3
	real(8) ,parameter :: sigmax =100.d0
	real(8) :: sig(0:L)

		do i=1,10
		sig(i) = sigmax* ((dble(i)-1.0d0/2.0d0)/L)**m
	
	write(*,*) sigmax,sig(i) 
 		enddo
write(*,*) 1/3, dble(1/3), dble(1)/dble(3),1.0d0/3.0d0, (10.0d0-1.0d0/2.0d0)/10.0d0,3.0d0/7.0d0
		endprogram
