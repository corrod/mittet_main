! ficticious E'(t') to diffusive frequency domain E(ω), using DFT
! frequency green function GX_w(ω)
! DFT
program f_to_d!(ns)
	use const_para
	implicit none

	integer :: nd,ios
	real(8), allocatable :: inp1_r(:),inp1_i(:),inp2_r(:),inp2_i(:),t1(:),t2(:)
	real(8) :: w(:) !窓関数

	integer :: n , l
	integer,intent(in) :: ns !sampling数
	complex(kind(0d0)),allocatable ::EX_w(:) !EX_w(0:ns-1) !周波数領域のEX
	complex(kind(0d0)),allocatable ::EX_f(:) !EX_f(0:ns-1) !ficticiousのE'x
	complex(kind(0d0)),allocatable ::JX_w(:) !JX_w(0:ns-1) !ficticiousのJ'x
	complex(kind(0d0)),allocatable ::JX_f(:) !JX_f(0:ns-1)
	complex(kind(0d0)),allocatable ::GX_w(:) !GX_w(0:ns-1) !diffusive domain Green's function
	character(3) :: name
!	omega=2.0d0*pi*k/ns


 !!! データの読み込み------------------------------------------------------------

    !EXファイル（データ）の長さNDを調べる-------------------------------------
    open(51,file='inp1.dat',action='read')
      nd=0
    do
        read(51,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(51)

     !配列確保
    allocate( w(0:nd-1),t1(0:nd-1),inp1_r(0:nd-1),inp1_i(0:nd-1),t2(0:nd-1),inp2_r(0:nd-1),inp2_i(0:nd-1),EX_w(0:nd-1),EX_f(0:nd-1),JX_w(0:nd-1),JX_f(0:nd-1),GX_w(0:nd-1) )

     !DFTするEXデータの読み込み
    open(51,file='inp1.dat',action='read')
      do i=1,nd
      read(51,*) t1(i), inp1_r(i), inp2_i(i)
      enddo
    close(51)

    EX_f(0:ns-1) = inp1_r(0:ns-1) + (0.0d0,1.0d0)*inp1_i(0:ns-1)


    !窓関数をかける----------------
    call window_cos(nd,w)

	!taper かけて
    do i=0,nd-1
    	Ex_f(i) = Ex_f(i) * w(i)
	enddo

	!JXファイル（データ）の長さNDを調べる------------------------------------
    open(51,file='inp2.dat',action='read')
      nd=0
    do
        read(51,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(51)


     !DFTするJXデータの読み込み
    open(51,file='inp2.dat',action='read')
      do i=1,nd
      read(51,*) t2(i), inp2_r(i), inp2_i(i)
      enddo
    close(51)


    JX_f(0:nd-1) = inp2_r(0:nd-1) + (0.0d0,1.0d0)*inp2_i(0:nd-1)


    !窓関数をかける cos_window-----------------
    call window_cos(nd,w)

	!taper かけて
    do i=0,nd-1
    	Jx_f(i) = Jx_f(i) * w(i)
	enddo


	EX_w(0:ns-1) = 0.0d0
	JX_w(0:ns-1) = 0.0d0

!DFT開始----------------------------------------------------------------------------
    ns = nd  !サンプリング数


	do k=0,ns-1  !周波数用ループ

		do n=0,ns-1 !時間用ループ

		EX_w(k) = EX_w(k) &
				+ EX_f(n) * exp( sqrt(2.0d0*pi*omega0*k/dble(ns)) * (I_u-1.0d0) * n )

		JX_w(k) = JX_w(k) &
				+ exp( -2.0d0*omega0/I_u/(2.0d0*pi*k/dble(ns)) ) &
				* JX_f(n) * exp( sqrt(2.0d0*pi*omega0*k/dble(ns)) * (I_u-1.0d0) * n )

		enddo

! 		JX_w(0) = 2.0d0 * omega0  !!!要確認

		GX_w(ｋ) = EX_w(ｋ) / JX_w(k)  !JX_w /= 0

	enddo

!出力-------------------------------------------------------------------------------------
! 		write(name,'(I3)') l  受信点位置とかがいいかも
		!ある点での周波数領域EX_w
! 		open(50,file='EX_w'//name/'.d')
		open(50,file='out1.dat')

		do k=0,ns-1
			write(50,*) k * 2.0d0*pi/ns, real(EX_w(k)),aimag(EX_w(k))
		enddo
		close(50)

		!ある点での周波数領域JX_w
! 		open(51,file='JX_w'//name/'.d')
		open(51,file='out2.dat')

		do k=0,ns-1
			write(51,*) k * 2.0d0*pi/ns, real(JX_w(k)),aimag(JX_w(k))
		enddo
		close(51)


		!ある点での周波数領域グリーン関数
! 		open(52,file='GX_w'//name/'.d')
		open(52,file='out3.dat')

		do k=0,ns-1
			write(52,*) k * 2.0d0*pi/ns, real(GX_w(k)),aimag(GX_w(k))
		enddo
		close(52)


deallocate ( w,t1,t2,inp1_r,inp1_i,inp2_r,inp2_i,EX_w,EX_f,JX_w,JX_f,GX_w )
end program f_to_d










! !!by k
! subroutine f_to_d_matrix
! 	use const_para
! 	implicit none
! 		integer :: s !sampling number 2**○
! 		do j=1,s
! 			do k=1,s
! 				A(j,k) = exp(-(2.0d0*pi*sqrt((j-1)*s*t)*(k-1)/dble(s) ) )   *exp(I_u*(2.0d0*pi*sqrt((j-1)*s*t)*(k-1)/dble(s) ))
! 			enddo
! 		enddo
! end subroutine f_to_d_matrix

! !by i
! subroutine laplace_fft
! 	use const_para
! 	implicit none
! 		integer :: n
! 		integer :: it
! ! 		integer :: istep!!
! 		real(8) :: om
! 		complex(kind(0d0)) :: Ex_w(nstep)
! 		complex(kind(0d0)) :: Jx_w(nstep)
! 		complex(kind(0d0)) :: Gx_w(nstep)

! 		om =2.0d0*pi/it/dt
! ! 		om =2.0d0*pi/nstep/dt
! ! 		om =2.0d0*pi/istep/dt
! 		t0=pi/fmax_w
! 		beta=pi*fmax**2.0d0


! 	do n=1,it  !what's it? !0~? 1~?

! 		do k=1,it
! 			EX_w(n) = EX_w(n) &
! 					 + EX_f(k)*dt *exp(-sqrt(omega0*om*n)*k*dt) *exp(I_u*sqrt(omega0*om*n)*k*dt)

! 			JX_w(n) = JX_w(n) &
! 					+ sqrt(-2.0d0*omega0/I_u/om/dble(n)) * JX_f(k)*dt *exp(-sqrt(omega0*om*n)*k*dt) *exp(I_u*sqrt(omega0*om*n)*k*dt)
! 		enddo
! 			JX_w(0) = 2.0d0 * omega0  !!!要確認

! 			GX_w(n) = EX_w(n) / JX_w(n)

! 	enddo
! end subroutine laplace_fft



! subroutine convolution_GJ_to_E
! 	use const_para
! 	implicit none
! 		complex(kind(0d0)) :: in_G
! 		complex(kind(0d0)) :: in_J
! 		complex(kind(0d0)) :: in_EF
! 		complex(kind(0d0)) :: out_G
! 		complex(kind(0d0)) :: out_J
! 		complex(kind(0d0)) :: out_ET

! end subroutine convolution_GJ_to_E














