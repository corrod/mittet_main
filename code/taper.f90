!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FFT用 窓関数taper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program taper
	use const_para
		implicit none
		integer, parameter :: m=5 !taper幅
		integer :: n,nd,ios
		real(8), allocatable :: inp1(:)
		real(8), allocatable :: temp1(:)
		real(8), allocatable :: w(:) !window function重み


		!データ長計算--------------------------------------------
		open(51,file='inp1.dat')
	    nd=0
   		do
	        read(51,'(f12.0)',iostat=ios)
	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
	        nd=nd+1
	    enddo
	    close(51)
	    allocate(inp(1:nd), temp(1:nd), w(:))
	

	    !cos window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    do i=1,nd
	    	w(i) = cos(acos(-1.0d0)*0.5d0*dble(i-20)/12.0d0)
		    write(*,*) w(i)
	    enddo


	    !Bartlett window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    do i=1,m/2
		    w(i) = (m/2-i)*2,0d0/m !n<=m/2
	    enddo

	    do i=m/2+1,nd
		    w(i) = 0.0d0
	    enddo


	    !Hanning window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    do i=1,m/2
		    w(i) = cos(pai*i/m)      !n<=m/2
	    enddo

	    do i=m/2+1,nd
		    w(i) = 0.0d0
	    enddo


	    !Hamming window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    do i=1,m/2
		    w(i) = 0.54d0 + 0.46d0*cos(2.0d0*pai*i/m)      !n<=m/2
	    enddo

	    do i=m/2+1,nd
		    w(i) = 0.0d0
	    enddo


		!Blackman window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		do i=1,m/2
		    w(i) = 0.42d0 + 0.5d0*cos(2.0d0*pai*i/m) - 0.08d0*cos(4.0d0*pai*i/m)      !n<=m/2
	    enddo

	    do i=m/2+1,nd
		    w(i) = 0.0d0
	    enddo


	    !taperかけるデータの読み込み----------------------
	    open(51,file='inp1.dat')
	    do i=1,nd
		    read(51,*) t(i), inp(i)
	    enddo
	    close(51)
		
		!taper かけて
	    do i=1,nd
	    	temp1(i) = inp1(i) * w(i)
    	enddo
    	!出力
    	open(51,file='inp1.dat')
    	do i=1,nd
    		write(51,*) t(i), temp1(i) 
		enddo

	end program taper