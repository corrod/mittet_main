!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FFT用 窓関数taper
! 範囲要確認
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!cos window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_cos(nd,w)
		use const_para
		implicit none
		integer, intent(in)  :: nd !taper幅
		real(8) :: w(0:nd-1) !window function重み

	    do i=0,nd-1
	    	w(i) = cos(acos(-1.0d0)*0.5d0*dble(i-20)/12.0d0)
! 		    write(*,*) w(i)
	    enddo
end subroutine window_cos

!Bartlett window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_bartlett(nd,w)
		use const_para
		implicit none
		integer, parameter   :: m=5 !taper幅
		integer, intent(in)  :: nd !taper幅
		real(8) :: w(0:nd-1) !window function重み

	    do i=1,m/2
		    w(i) = (dble(m)/2.0d0-dble(i))*2.0d0/dble(m) !n<=m/2
	    enddo
end subroutine window_bartlett

!Hanning window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_hanning(nd,w)
	use const_para
		implicit none
		integer, parameter   :: m=5 !taper幅
		integer, intent(in)  :: nd !taper幅
		real(8) :: w(0:nd-1) !window function重み

	    do i=1,m/2
		    w(i) = cos(pi*dble(i)/dble(m) )     !n<=m/2
	    enddo
end subroutine window_hanning

!Hamming window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_hamming(nd,w)
	use const_para
		implicit none
		integer, parameter   :: m=5 !taper幅
		integer, intent(in)  :: nd !taper幅
		real(8) :: w(0:nd-1) !window function重み

	    do i=1,m/2
		    w(i) = 0.54d0 + 0.46d0*cos(2.0d0*pi*dble(i)/dble(m))      !n<=m/2
	    enddo
end subroutine window_hamming

!Blackman window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine window_blackman(nd,w)
	use const_para
		implicit none
		integer, parameter   :: m=5 !taper幅
		integer, intent(in)  :: nd !taper幅
		real(8) :: w(0:nd-1) !window function重み

		do i=1,m/2
		    w(i) = 0.42d0 + 0.5d0*cos(2.0d0*pi*dble(i)/dble(m)) - 0.08d0*cos(4.0d0*pi*dble(i)/dble(m))      !n<=m/2
	    enddo
end subroutine window_blackman














! subroutine cos_window
! 	use const_para
! 		implicit none
! 		integer, parameter   :: m=5 !taper幅
! 		integer              :: n,nd,ios
! 		real(8), allocatable :: inp1(:)
! 		real(8), allocatable :: temp1(:)
! 		real(8), allocatable :: w(:) !window function重み

! !データ長計算--------------------------------------------
! 		open(51,file='inp1.dat')
! 	    nd=0
!    		do
! 	        read(51,'(f12.0)',iostat=ios)
! 	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! 	        nd=nd+1
! 	    enddo
! 	    close(51)
! 	    allocate(inp1(1:nd), temp1(1:nd), w(:))

! 	    !cos window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	    do i=1,nd
! 	    	w(i) = cos(acos(-1.0d0)*0.5d0*dble(i-20)/12.0d0)
! 		    write(*,*) w(i)
! 	    enddo
!  !taperかけるデータの読み込み----------------------
! 	    open(51,file='inp1.dat')
! 	    do i=1,nd
! 		    read(51,*) t(i), inp(i)
! 	    enddo
! 	    close(51)

! 		!taper かけて
! 	    do i=1,nd
! 	    	temp1(i) = inp1(i) * w(i)
!     	enddo
!     	!出力
!     	open(51,file='inp1.dat')
!     	do i=1,nd
!     		write(51,*) t(i), temp1(i)
! 		enddo

! 	end subroutine cos_window





! subroutine bartlett_window
! 	use const_para
! 		implicit none
! 		integer, parameter   :: m=5 !taper幅
! 		integer              :: n,nd,ios
! 		real(8), allocatable :: inp1(:)
! 		real(8), allocatable :: temp1(:)
! 		real(8), allocatable :: w(:) !window function重み

! !データ長計算--------------------------------------------
! 		open(51,file='inp1.dat')
! 	    nd=0
!    		do
! 	        read(51,'(f12.0)',iostat=ios)
! 	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! 	        nd=nd+1
! 	    enddo
! 	    close(51)
! 	    allocate(inp(1:nd), temp(1:nd), w(:))
! 	    !Bartlett window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	    do i=1,m/2
! 		    w(i) = (dble(m)/2.0d0-dble(i))*2,0d0/dble(m) !n<=m/2
! 	    enddo

! 	    do i=m/2+1,nd
! 		    w(i) = 0.0d0
! 	    enddo
!  !taperかけるデータの読み込み----------------------
! 	    open(51,file='inp1.dat')
! 	    do i=1,nd
! 		    read(51,*) t(i), inp(i)
! 	    enddo
! 	    close(51)

! 		!taper かけて
! 	    do i=1,nd
! 	    	temp1(i) = inp1(i) * w(i)
!     	enddo
!     	!出力
!     	open(51,file='inp1.dat')
!     	do i=1,nd
!     		write(51,*) t(i), temp1(i)
! 		enddo

! 	end subroutine bartlett_window





! subroutine hanning_window
! 	use const_para
! 		implicit none
! 		integer, parameter   :: m=5 !taper幅
! 		integer              :: n,nd,ios
! 		real(8), allocatable :: inp1(:)
! 		real(8), allocatable :: temp1(:)
! 		real(8), allocatable :: w(:) !window function重み

! !データ長計算--------------------------------------------
! 		open(51,file='inp1.dat')
! 	    nd=0
!    		do
! 	        read(51,'(f12.0)',iostat=ios)
! 	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! 	        nd=nd+1
! 	    enddo
! 	    close(51)
! 	    allocate(inp(1:nd), temp(1:nd), w(:))
! 	    !Hanning window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	    do i=1,m/2
! 		    w(i) = cos(pi*dble(i)/dble(m) )     !n<=m/2
! 	    enddo

! 	    do i=m/2+1,nd
! 		    w(i) = 0.0d0
! 	    enddo
!  !taperかけるデータの読み込み----------------------
! 	    open(51,file='inp1.dat')
! 	    do i=1,nd
! 		    read(51,*) t(i), inp(i)
! 	    enddo
! 	    close(51)

! 		!taper かけて
! 	    do i=1,nd
! 	    	temp1(i) = inp1(i) * w(i)
!     	enddo
!     	!出力
!     	open(51,file='inp1.dat')
!     	do i=1,nd
!     		write(51,*) t(i), temp1(i)
! 		enddo

! 	end subroutine hanning_window





! subroutine hamming_window
! 	use const_para
! 		implicit none
! 		integer, parameter   :: m=5 !taper幅
! 		integer              :: n,nd,ios
! 		real(8), allocatable :: inp1(:)
! 		real(8), allocatable :: temp1(:)
! 		real(8), allocatable :: w(:) !window function重み
! !データ長計算--------------------------------------------
! 		open(51,file='inp1.dat')
! 	    nd=0
!    		do
! 	        read(51,'(f12.0)',iostat=ios)
! 	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! 	        nd=nd+1
! 	    enddo
! 	    close(51)
! 	    allocate(inp(1:nd), temp(1:nd), w(:))
! 	    !Hamming window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 	    do i=1,m/2
! 		    w(i) = 0.54d0 + 0.46d0*cos(2.0d0*pi*dble(i)/dble(m))      !n<=m/2
! 	    enddo

! 	    do i=m/2+1,nd
! 		    w(i) = 0.0d0
! 	    enddo
! !taperかけるデータの読み込み----------------------
! 	    open(51,file='inp1.dat')
! 	    do i=1,nd
! 		    read(51,*) t(i), inp(i)
! 	    enddo
! 	    close(51)

! 		!taper かけて
! 	    do i=1,nd
! 	    	temp1(i) = inp1(i) * w(i)
!     	enddo
!     	!出力
!     	open(51,file='inp1.dat')
!     	do i=1,nd
!     		write(51,*) t(i), temp1(i)
! 		enddo

! 	end subroutine hamming_window






! subroutine blackman_window
! 	use const_para
! 		implicit none
! 		integer, parameter   :: m=5 !taper幅
! 		integer              :: n,nd,ios
! 		real(8), allocatable :: inp1(:)
! 		real(8), allocatable :: temp1(:)
! 		real(8), allocatable :: w(:) !window function重み


! !データ長計算--------------------------------------------
! 		open(51,file='inp1.dat')
! 	    nd=0
!    		do
! 	        read(51,'(f12.0)',iostat=ios)
! 	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! 	        nd=nd+1
! 	    enddo
! 	    close(51)
! 	    allocate(inp(1:nd), temp(1:nd), w(:))


! 		!Blackman window!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 		do i=1,m/2
! 		    w(i) = 0.42d0 + 0.5d0*cos(2.0d0*pi*dble(i)/dble(m)) - 0.08d0*cos(4.0d0*pi*dble(i)/dble(m))      !n<=m/2
! 	    enddo

! 	    do i=m/2+1,nd
! 		    w(i) = 0.0d0
! 	    enddo


! !taperかけるデータの読み込み----------------------
! 	    open(51,file='inp1.dat')
! 	    do i=1,nd
! 		    read(51,*) t(i), inp(i)
! 	    enddo
! 	    close(51)

! 		!taper かけて
! 	    do i=1,nd
! 	    	temp1(i) = inp1(i) * w(i)
!     	enddo
!     	!出力
!     	open(51,file='inp1.dat')
!     	do i=1,nd
!     		write(51,*) t(i), temp1(i)
! 		enddo

! 	end subroutine blackman_window