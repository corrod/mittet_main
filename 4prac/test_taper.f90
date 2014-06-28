! program taper
! 		implicit none
! 		integer,parameter  :: nd=1000 !taper幅
! 		real(8) :: w(0:nd-1) !window function重み
! 		integer ::i
! 		open(1,file='taper.d')
! 	    do i=0,nd-1
! 	    	w(i) = cos(acos(-1.0d0)*0.5d0*dble(i-10)/12.0d0)
! 		    write(1,*) i,w(i)
! 	    enddo
! 	    close(1)
! 	end program taper

! program hanning
! ! 	use const_para
! 		implicit none
! 		integer i
! ! 		integer, parameter   :: m=1000 !taper幅
! 		integer, parameter  :: nd=10000 !taper幅
! 		real(8) :: w(0:nd-1) !window function重み
! 		open(1,file='taper.d')

! ! 	    do i=1,nd!m/2
! ! 		    w(i) = cos(3.1415d0*dble(i)/dble(m) )
! ! 		    write(1,*) i,w(i)

! 	    do i=1,nd
! 	    	w(i) = 1-(cos(3.1415d0*i/nd))**2
! ! 			w(i) = 0.5d0 + 0.5d0*cos(2.0d0*3.1415d0*i/nd)
! 		    write(1,*) i,w(i)
! 		     !n<=m/2
! 	    enddo
! 	    close(1)

! end program hanning

! program hamming
! ! 	use const_para
! 		implicit none
! 		integer i
! ! 		integer, parameter   :: m=1000 !taper幅
! 		integer, parameter  :: nd=10000 !taper幅
! 		real(8) :: w(0:nd-1) !window function重み
! 		open(1,file='taper.d')

! ! 	    do i=1,nd!m/2
! ! 		    w(i) = cos(3.1415d0*dble(i)/dble(m) )
! ! 		    write(1,*) i,w(i)

! 	    do i=1,nd
! 	    	w(i) = 0.54d0-0.46d0*cos(2.0d0*3.14d0*i/nd)
! ! 			w(i) = 0.5d0 + 0.5d0*cos(2.0d0*3.1415d0*i/nd)
! 		    write(1,*) i,w(i)
! 		     !n<=m/2
! 	    enddo
! 	    close(1)

! end program hamming

program flattop
! 	use const_para
		implicit none
		integer i
		integer, parameter  :: nd=10000 !taper幅
		real(8) :: w1(0:nd-1) !window function重み
		real(8) :: w(0:nd-1)
		open(1,file='taper.d')
	    do i=0,nd-1
! 	    	w(i) = 1-(cos(3.14*dble(i)/dble(nd-1)))**2.0d0
	    	w(i) = 0.54d0-0.46d0*cos(2.0d0*3.14d0*dble(i)/dble(nd-1))!hamming
! 	    	w(i) = 0.54d0-0.46d0*cos(2.0d0*3.14d0*dble(i)/dble(nd-1))

! 	    	w(i) = cos(acos(-1.0d0)*0.5d0*dble(i-20)/12.0d0)
! 	    	    w(i) = 0.5d0 - 0.42d0*cos(2.0d0*3.141592d0*dble(i)/dble(nd-1)) - 0.08d0*cos(4.0d0*3.141592d0*dble(i)/dble(nd-1))
! 	    	w(i) = 0.54d0-0.46d0*cos(2.0d0*3.14d0*dble(i)/dble(nd-1))
! 	    	w1(i) = (4.0d0*3.14d0*(dble(i)-dble(nd-1)/2.0d0)/dble(nd-1.0d0))
! 	    	w(i) = (0.54d0+0.46d0*cos(0.5d0*w1(i)))*sin(w1(i))/w1(i)
 write(1,*) i,w(i)
	    enddo
	    close(1)

end program flattop