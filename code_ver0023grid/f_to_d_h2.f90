!///////////////////////////////////////////////////////////////////////////////
!conjugate version     H
!////////////////////////////////////////////////////////////////////////////
! ficticious H'(t') to diffusive frequency domain H(ω), using DFT, FFT
! frequency green function GXh_w(ω)
!
!JZ_w GXh_w ひとつめNAN  >> JZ(0) =2omega_0 　
!Je(istep) = dt*etaxx(x0,y0,z0)*signal(istep) /dx/dy/dz
!JZ_f = Jh(istep) = signal(istep)*dt / myu(x0,y0,z0) /dx/dy/dz☓
!JZ_f = signal(istep) ◎
!
!fwi3d_cpm_function 行1430~参照
!//////////////////////////////////////////////////////////////////////////

program f_to_d_h
	use const_para
	implicit none

	integer :: nd,ios
	real(8), allocatable :: inp1_r(:),inp1_i(:),inp2_r(:),inp2_i(:),t1(:),t2(:)
	real(8), allocatable :: w(:) !窓関数
	real(8) :: om

	integer :: n !, l
	complex(kind(0d0)),allocatable ::Hz_w(:) !Hz_w(0:nd-1) !周波数領域のHz
	complex(kind(0d0)),allocatable ::Hz_f(:) !Hz_f(0:nd-1) !ficticiousのH'z
	complex(kind(0d0)),allocatable ::JZ_w(:) !JZ_w(0:nd-1) !diffusiveのJz
	complex(kind(0d0)),allocatable ::JZ_f(:) !JZ_f(0:nd-1)
	complex(kind(0d0)),allocatable ::GXh_w(:) !GXh_w(0:nd-1) !diffusive domain Green's function
    complex(kind(0d0)),allocatable :: inv_JZ_w(:)
	character(3) :: name
	!IDFT, IFFT用
	complex(kind(0d0)),allocatable :: Hz_t(:), JZ_t(:),GXh_t(:)
	complex(kind(0d0)),allocatable :: in1(:), in2(:), in3(:) !IFFT用
	complex(kind(0d0)),allocatable :: out1(:), out2(:), out3(:) !IFFT用
	integer(8) :: plan1, plan2, plan3

include 'fftw3.f'



!/////////////////////////////////////////////////////////////////////////////
! データの読み込み Hz_f
!/////////////////////////////////////////////////////////////////////////////
    !Hzファイル（データ）の長さNDを調べる-------------------------------------
    open(51,file='inp1.dat',action='read')
      nd=0
	    do
	        read(51,'(f12.0)',iostat=ios)
	        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
	         nd=nd+1
	    enddo
    close(51)


!///////////////////////////////////////////////////////////////////////////////
! 配列確保
!///////////////////////////////////////////////////////////////////////////////
    allocate(t1(0:nd-1),inp1_r(0:nd-1),inp1_i(0:nd-1),t2(0:nd-1),inp2_r(0:nd-1),inp2_i(0:nd-1))
    allocate(w(0:nd-1),Hz_w(0:nd-1),Hz_f(0:nd-1),JZ_w(0:nd-1),JZ_f(0:nd-1),GXh_w(0:nd-1),inv_JZ_w(0:nd-1))



!DFTするHzデータの読み込み
open(51,file='inp1.dat',action='read')
  do i=0,nd-1
    read(51,*) t1(i), inp1_r(i), inp1_i(i)
  Hz_f(i) = inp1_r(i) + (0.0d0,1.0d0)*inp1_i(i)
  enddo
close(51)


!////////////////////////////////////////////////////////////////////////////
!Hz_fに窓関数をかける hamming window
!///////////////////////////////////////////////////////////////////////////
!     call window_hamming(nd,w) !hamming 両端が0にはならない窓
! !     call window_hanning(nd,w) !hanning 両端が0になる窓
!       !taper かけて
!       do i=0,nd-1
!           write(8,*) i*dt,real(Hz_f(i)),aimag(Hz_f(i))!かける前
!           Hz_f(i) = Hz_f(i) * w(i)
!           write(9,*) i*dt,real(Hz_f(i)),aimag(Hz_f(i))!かけた後
!       enddo


!/////////////////////////////////////////////////////////////////////////////
! データの読み込み JZ_f
!/////////////////////////////////////////////////////////////////////////////
    !JZファイル（データ）の長さNDを調べる------------------------------------
    open(51,file='inp2.dat',action='read')
      nd=0
        do
            read(51,'(f12.0)',iostat=ios)
            if (ios<0) exit !ファイルの末尾にきたらループを抜ける
             nd=nd+1
        enddo
    close(51)

     !DFTするJZデータの読み込み
    open(51,file='inp2.dat',action='read')
      do i=0,nd-1
          read(51,*) t2(i), inp2_r(i), inp2_i(i)
      JZ_f(i) = inp2_r(i) + (0.0d0,1.0d0)*inp2_i(i)
      enddo
    close(51)


!////////////////////////////////////////////////////////////////////////////
!JZ_fに窓関数をかける hamming window
! !///////////////////////////////////////////////////////////////////////////
!     call window_hamming(nd,w) !hamming 両端が0にはならない窓
! !     call window_hanning(nd,w) !hanning 両端が0になる窓
!       do i=0,nd-1
!       write(7,*) w(i)
!       enddo
!   !taper かけて
!       do i=0,nd-1
!           write(10,*) i*dt,real(JZ_f(i)),aimag(JZ_f(i))!かける前出力
!           JZ_f(i) = JZ_f(i) * w(i)
!           write(11,*) i*dt,real(JZ_f(i)),aimag(JZ_f(i))!かけた後出力
!       enddo




!/////////////////////////////////////////////////////////////////////////////////
! DFT開始 ficticious to diffusive freq
!/////////////////////////////////////////////////////////////////////////////////
    write(*,*) '*********************        DFT start       ********************'


    om   = 2.d0*pi/dble(nd)/dt

    Hz_w(0:nd-1) = 0.0d0
    JZ_w(0:nd-1) = 0.0d0
    GXh_w(0:nd-1) = 0.0d0

    do k=0,nd-1  !周波数用ループ
    Hz_w(k) = 0.0d0
    JZ_w(k) = 0.0d0
        do n=0,nd-1 !時間用ループ

        ! mittet(11)の係数参照
        Hz_w(k) = Hz_w(k) &
                + sqrt( - 2.0d0*omega0/I_u/om/k ) * Hz_f(n) * dt &
                * exp( (I_u - 1.0d0) * sqrt(omega0*om*k) *  n*dt )

        ! (11) from mittet J(x,omega) = sqrt(-2*omega0/i*omega)*J'(x,omega)
!         JZ_w(k) = JZ_w(k) &
!               + sqrt( -2.0d0*omega0/I_u/om/k ) * JZ_f(n) *dt &
!               * exp( sqrt(omega0*om*k) * (I_u-1.0d0) * n*dt)

        ! (11) from mittet  K(x,omega) = K'(x,omega)　　　
        JZ_w(k) = JZ_w(k) &
                + JZ_f(n) * dt &
                * exp( (I_u - 1.0d0) * sqrt(omega0*om*k) *  n*dt )

        enddo !n loop

!         Hz_w(0) = 2.0d0 * omega0  !!! 　　　

        !(C-11)
!       JZ_w(0) = 2.0d0 * omega0  !!! 　　　

    GXh_w(k) = Hz_w(k) / JZ_w(k)  !JZ_w /= 0

    enddo !k loop

!///////////////////////////////////////////////////////////////////////////////
! output
!//////////////////////////////////////////////////////////////////////////////
!       write(name,'(I3)') l  受信点位置とかがいいかも
        !ある点での周波数領域Hz_w
!       open(50,file='Hz_w'//name/'.d')
        open(60,file='out1.dat')
        do k=0,nd-1
            write(60,*) k*om/2.0d0/pi, real(Hz_w(k)),aimag(Hz_w(k))   !!!横軸周波数の書き方違うかも
        enddo
        close(60)

        !ある点での周波数領域JZ_w
!       open(61,file='JZ_w'//name/'.d')
        open(61,file='out2.dat')
        do k=0,nd-1
            write(61,*) k*om/2.0d0/pi, real(JZ_w(k)),aimag(JZ_w(k))
        enddo
        close(61)


        !ある点での周波数領域グリーン関数
!       open(62,file='GXh_w'//name/'.d')
        open(62,file='out3.dat')
        do k=0,nd-1
            write(62,*) k*om/2.0d0/pi, real(GXh_w(k)),aimag(GXh_w(k))
        enddo
        close(62)

!         open(63,file='inv_JZ_w.dat')
!         do k=0,nd-1
!             write(63,*) k*om, real(inv_JZ_w(k)),aimag(inv_JZ_w(k))
!         enddo
!         close(63)

        !絶対値__________________________
        open(70,file='out4.dat')
        do k=0,nd-1
            write(70,*) k*om/2.0d0/pi, abs(HZ_w(k))  !!!横軸周波数の書き方違うかも
        enddo
        close(70)

        !ある点での周波数領域JZ_w
!       open(61,file='JZ_w'//name/'.d')
        open(71,file='out5.dat')
        do k=0,nd-1
            write(71,*) k*om/2.0d0/pi, abs(JZ_w(k))
        enddo
        close(71)


        !ある点での周波数領域グリーン関数
!       open(62,file='GXh_w'//name/'.d')
        open(72,file='out6.dat')
        do k=0,nd-1
            write(72,*) k*om/2.0d0/pi, abs(GXh_w(K))
        enddo
        close(72)

!//////////////////////////////////////////////////////////////////////////////////
!
! IFFT    Frequency to time trandformation   JZ_w,Hz_w,GXh_w to JZ_t,Hz_t,GXh_w
!
!/////////////////////////////////////////////////////////////////////////////////
    write(*,*) '********************        IFFT start       ********************'

nd = (nd-1) * 2
! nd = nd * 2

    allocate( in1(0:nd-1), in2(0:nd-1), in3(0:nd-1) )
    allocate( Hz_t(1:nd),JZ_t(1:nd),GXh_t(1:nd) )
    allocate( out1(1:nd), out2(1:nd), out3(1:nd) )


    in1(0:nd-1) = 0.0d0
    in2(0:nd-1) = 0.0d0
    in3(0:nd-1) = 0.0d0
    Hz_t(1:nd) = 0.0d0
    JZ_t(1:nd) = 0.0d0
    GXh_t(1:nd) = 0.0d0
    out1(1:nd) = 0.0d0
    out2(1:nd) = 0.0d0
	out3(1:nd) = 0.0d0

!////////////////////////////////////////////////////////////////////////////
! HZ_w,Jz_w,GXh_w に窓関数をかける hamming window
!///////////////////////////////////////////////////////////////////////////
! !     call window_hamming(nd,w) !hamming 両端が0にはならない窓
!     call window_hanning(nd,w) !hanning 両端が0になる窓
!       do i=0,nd-1
!       write(7,*) w(i)
!       enddo
!   !taper かけて
!       do k=0,nd-1
!           write(10,*) k*om/2.0d0/pi,real(HZ_w(k)),aimag(HZ_w(k))!かける前出力
!              HZ_w(k) = HZ_w(k) * w(k)
!              JZ_w(k) = JZ_w(k) * w(k)
!              GXh_w(k) = GXe_w(k) * w(k)
!           write(11,*) k*om/2.0d0/pi,real(HZ_w(k)),aimag(HZ_w(k))!かけた後出力
!       enddo

	do k=0,nd/2
		in1(k) = Hz_w(k)
		in2(k) = JZ_w(k)
		in3(k) = GXh_w(k)
	enddo

    in1(nd/2+1:nd-1) = conjg(Hz_w(nd/2-1:1:-1))
    in2(nd/2+1:nd-1) = conjg(JZ_w(nd/2-1:1:-1))
    in3(nd/2+1:nd-1) = conjg(GXh_w(nd/2-1:1:-1))

    open(101,file='conjg_hzw.dat')
    do i=0,nd-1
    write(101,*) i, real(in1(i)), aimag(in1(i))
    enddo
    close(101)

    open(102,file='conjg_jzw.dat')
    do i=0,nd-1
    write(102,*) i, real(in2(i)), aimag(in2(i))
    enddo

    open(103,file='conjg_gxhw.dat')
    do i=0,nd-1
    write(103,*) i, real(in3(i)), aimag(in3(i))
    enddo

!////////////////////////////////////////////////////////////////////////////
! in1,in2,in3 に窓関数をかける hamming window
!///////////////////////////////////////////////////////////////////////////
! !     call window_hamming(nd,w) !hamming 両端が0にはならない窓
!     call window_hanning(nd,w) !hanning 両端が0になる窓
!       do i=0,nd-1
!       write(7,*) w(i)
!       enddo
!   !taper かけて
!       do k=0,nd-1
!              in1(k) = in1(k) * w(k)
!              in2(k) = in2(k) * w(k)
!              in3(k) = in3(k) * w(k)
!       enddo

!////////////////////////////////////////////////////////////
! make pland
!       FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
!////////////////////////////////////////////////////////////
	call dfftw_plan_dft_1d(plan1,nd,in1,out1,FFTW_BACKWARD,FFTW_ESTIMATE) !complex array入力
	call dfftw_plan_dft_1d(plan2,nd,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE)
	call dfftw_plan_dft_1d(plan3,nd,in3,out3,FFTW_BACKWARD,FFTW_ESTIMATE)

!     call dfftw_plan_dft_1d(plan1,nd,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE) !complex array入力
!     call dfftw_plan_dft_1d(plan2,nd,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE)
!     call dfftw_plan_dft_1d(plan3,nd,in3,out3,FFTW_FORWARD,FFTW_ESTIMATE)

!///////////////////////////////////////////////////////////
! carry out fourier trandformation
!///////////////////////////////////////////////////////////
	call dfftw_execute_dft(plan1,in1,out1) !　_dft
	call dfftw_execute_dft(plan2,in2,out2)
	call dfftw_execute_dft(plan3,in3,out3)


!///////////////////////////////////////////////////////////
! destroy plan
!///////////////////////////////////////////////////////////
	call dfftw_destroy_plan(plan1)
	call dfftw_destroy_plan(plan2)
	call dfftw_destroy_plan(plan3)


!////////////////////////////////////////////////////////////
! output
!////////////////////////////////////////////////////////////
	open(81,file='invGH.dat')
	open(82,file='invGJ.dat')
	open(83,file='invGG.dat')
    open(84,file='absHZ_t.dat')
    open(85,file='absJZ_t.dat')
    open(86,file='absGXh_t.dat')

	do n=1,nd
        !スケール /nd/dt*2.0d0 　　　
        out1(n) = out1(n)/nd/dt *2.0d0
        out2(n) = out2(n)/nd/dt *2.0d0
        out3(n) = out3(n)/nd/dt *2.0d0
         !スケール / nd 　　　
!         out1(n) = out1(n)/nd !E
!         out2(n) = out2(n)/nd !J
!         out3(n) = out3(n)/nd !GX_t

		write(81,*) n*dt, real(out1(n)), aimag(out1(n))
		write(82,*) n*dt, real(out2(n)), aimag(out2(n))
		write(83,*) n*dt, real(out3(n)), aimag(out3(n))
        write(84,*) n*dt, abs(out1(n))
        write(85,*) n*dt, abs(out2(n))
        write(86,*) n*dt, abs(out3(n))

        GXh_t(n) = out3(n)
	enddo
	close(81)
	close(82)
    close(83)
    close(84)
    close(85)
	close(86)

	deallocate( w,t1,t2,inp1_r,inp1_i,inp2_r,inp2_i,Hz_w,Hz_f,JZ_w,JZ_f,GXh_w,inv_JZ_w )
	deallocate( in1,in2,in3,out1,out2,out3,Hz_t,JZ_t,GXh_t )

end program f_to_d_h
















!///////////////////////////////////////////////////////////////////////////////////////
!
! IDFT     Frequency to time trandformation JZ_w,Hz_w,GXh_w to JZ_t,Hz_t,GXh_w
!
!///////////////////////////////////////////////////////////////////////////////////////
!   write(*,*) '*********************        IDFT start       ************t********'

!   Hz_t(0:nd-1) = 0.0d0
!   JZ_t(0:nd-1) = 0.0d0
!   GXh_w(0:nd-1) = 0.0d0

!   do k=0,nd-1
!       do n=0,nd-1
!       Hz_t(k) = Hz_t(k) &
!               + Hz_w(n) * exp(-I_u*2.0d0*pi*k*n/nd) /nd/dt *2.0d0
!       JZ_t(k) = JZ_t(k) &
!               + JZ_w(n) * exp(-I_u*2.0d0*pi*k*n/nd) /nd/dt *2.0d0
!       GXh_w(k) = GXh_w(k) &
!               + GXh_w(n) * exp(-I_u*2.0d0*pi*k*n/nd) /nd/dt *2.0d0
!       enddo
!   enddo

!       open(71,file='invGH.dat')
!       open(72,file='invGJ.dat')
!       open(73,file='invGG.dat')
!   do k=0,nd-1
!       write(71,*) k*dt, real(Hz_t(k)), aimag(Hz_t(k))
!       write(72,*) k*dt, real(JZ_t(k)), aimag(JZ_t(k))
!       write(73,*) k*dt, real(GXh_w(k)), aimag(GXh_w(k))
!   enddo
!       close(71)
!       close(72)
!       close(73)

!   deallocate( w,t1,t2,inp1_r,inp1_i,inp2_r,inp2_i,Hz_w,Hz_f,JZ_w,JZ_f,GXh_w )
!   deallocate( Hz_t, JZ_t, GXh_w )

! end program f_to_d













!!ficticiou wave domain からdiffusive frequency using FFT //////////////////////////
! FFTW
!入出力の列数に注意
!deallocateするとプログラム止まる問題ありls
!///////////////////////////////////////////////////////////////////////////////////

! program f_to_d
! 	use condt_para
! 		implicit none

!     		integer                         :: n,nd,ios
!     		!real(8)                        :: omega
!         real(8), allocatable            :: inp(:), t(:), f(:), p(:)
!         complex(kind(0d0)), allocatable :: c(:)
!     	!	complex(kind(0d0)), allocatable :: Gxn(:,:,:),Gyn(:,:,:),Gzn(:,:,:)
!     	!	real(8), allocatable :: jh(:)
!        ! character(8) :: inp1
!        ! character(8) :: out1,out2,out3
!        ! character(5), parameter :: inp2='jh.d'
!        ! character(12),parameter :: out4='jh_as.out',out5='jh_ps.out',out6='jh_fft.out'
!         integer :: plan(8)


! 	! FFTW3を呼び出すのに必要なヘッダーファイルを include する
!     include 'fftw3.f'

!     !!! 開始------------------------------------------------------------
!     !ファイル（データ）の長さNDを調べる
!     open(51,file='inp1.dat',action='read')
!       nd=0
!     do
!         read(51,'(f12.0)',iostat=ios)
!         if (ios<0) exit !ファイルの末尾にきたらループを抜ける
!          nd=nd+1
!     enddo
!     close(51)


!   	!データの長さをnの２乗になるように決めて
!     !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
!     n = 2**int(log(dble(nd))/ log(2.0d0)+0.5d0)

!     !配列確保
!     allocate(t(1:n),inp(1:n),c(1:n/2),f(1:n/2),p(1:n/2))


!     !fftするデータの読み込み
!     open(51,file='inp1.dat',action='read')
!       do i=1,nd
!       read(51,*) t(i), inp(i)
!       enddo
!     close(51)


!    	!不足分は0パッディング
!     inp(nd+1:n) = 0.0d0


!     !fftwの実行
!     call dfftw_plan_dft_r2c_1d(plan,n,inp,c,fftw_estimate)
!     call dfftw_execute(plan)
!     call dfftw_destroy_plan(plan)


!     ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
!     f= abs(c)
!     p= atan2(aimag(c),dble(c))


! 	! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
!     open(61,file='out1.dat')!,status='replace',action='write')
!       do i=1,n/2
!       write(61,*) dble(i-1)/(dble(n)*dt), f(i)
!       enddo
!     close(61)


!     open(61,file='out2.dat')!,status='replace',action='write')
!       do i=1,n/2
!       write(61,*) dble(i-1)/(dble(n)*dt), p(i)
!       enddo
!     close(61)


!      !フーリエ変換後のdiffusive freqency domainのinp
!     open(61,file='out3.dat')!,status='replace',action='write')
!       do i=1,n/2
!       write(61,*) dble(i-1)/(dble(n)*dt), real(c(i)), aimag(c(i))
!       enddo
!     close(61)

!     end program f_to_d

! 	!!!jh-------------------------------------------------------------------
! 	!ファイル（データ）の長さNDを調べる
! !    open(52,file=inp2,action='read')
! !    nd=0
! !    do
! !        read(52,'(f12.0)',iostat=ios)
! !        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
! !         nd=nd+1
! !!    enddo
!  !   close(52)

!  	!データの長さをnの２乗になるように決めて
!     !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
! !     n=2**int(log(dble(nd))/ log(2.0d0)+0.5d0)

!  !    allocate(t(1:n),jh(1:n),c(1:n/2),f(1:n/2),p(1:n/2))
! !    write(*,*) 'データ数nd,nの２乗数',nd,n

!  !   if (nd>n) then
!  !!   write(*,*) "nd must < n"
  !  endif

  !  allocate(jh(1:n))

    !fftするデータの読み込み
!	open(52,file=inp2,action='read')          !!!!!!!!!!!!!!  n>ndが崩れるときがある
 !   read(52,*) jh(1:nd)
  !  close(52)


 	!不足分は0パッディング
  !  jh(nd+1:n)=0.0d0

    !fftwの実行
 !   call dfftw_plan_dft_r2c_1d(plan,n,jh,c,fftw_estimate)
  !  call dfftw_execute(plan)
   ! call dfftw_destroy_plan(plan)

    ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
    !f= abs(c)
    !p= atan2(aimag(c),dble(c))

	! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
 !   open(61,file=out4,status='replace',action='write')
  !  do i=1,n/2
   ! write(61,*) f(i)
    !enddo
 !   close(61)

  !  open(61,file=out5,status='replace',action='write')
  !  do i=1,n/2
  !  write(61,*) p(i)
  !  enddo
  !  close(61)

     !フーリエ変換後のdiffusive freqency domainのjh
        !√-2omega0/iomeagaをかける
!    jh(1:n) = sqrt(-2.0d0*omega0/(i*omega)) * c(0:n/2)
 !   open(61,file=out6,status='replace',action='write')
  !  do i=1,n/2
   ! write(61,*) (i-1)/(n*5.0d-4),real(c(i)),aimag(c(i))
    !enddo
    !close(61)


!!!グリーン関数を求める-------------------------------------------------
!	Gxn(:,:,:) = Ex / jh
!	Gyn(:,:,:) = Ey / jh
!	Gzn(:,:,:) = Ez / jh












