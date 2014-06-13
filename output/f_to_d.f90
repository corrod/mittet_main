!!ficticiou wave domain からdiffusive frequency !!!!!!!
!
!入出力の列数に注意
!deallocateするとプログラム止まる問題ありls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program f_to_d
	use const_para
		implicit none

    		integer :: n,nd,ios
    		real(8) :: omega 
    	!	complex(kind(0d0)), allocatable :: Gxn(:,:,:),Gyn(:,:,:),Gzn(:,:,:)
    	!	real(8), allocatable :: jh(:)
        real(8), allocatable :: inp(:), t(:), f(:), p(:)
  	   	complex(kind(0d0)), allocatable :: c(:)
       ! character(8) :: inp1
       ! character(8) :: out1,out2,out3
       ! character(5), parameter :: inp2='jh.d'
       ! character(12),parameter :: out4='jh_as.out',out5='jh_ps.out',out6='jh_fft.out'
        integer :: plan(8)


	! FFTW3を呼び出すのに必要なヘッダーファイルを include する
    include 'fftw3.f'
 
    !!! 開始------------------------------------------------------------
    !ファイル（データ）の長さNDを調べる
    open(51,file='inp1.dat',action='read')
      nd=0
    do
        read(51,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(51)


  	!データの長さをnの２乗になるように決めて
    !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
    n = 2**int(log(dble(nd))/ log(2.0d0)+0.5d0)

    !配列確保
    allocate(t(1:n),inp(1:n),c(1:n/2),f(1:n/2),p(1:n/2))

	
    !fftするデータの読み込み
    open(51,file='inp1.dat',action='read')
      do i=1,nd
      read(51,*) t(i), inp(i)
      enddo
    close(51)


   	!不足分は0パッディング
    inp(nd+1:n) = 0.0d0


    !fftwの実行
    call dfftw_plan_dft_r2c_1d(plan,n,inp,c,fftw_estimate)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)


    ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
    f= abs(c)
    p= atan2(aimag(c),dble(c))


	! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
    open(61,file='out1.dat')!,status='replace',action='write')
      do i=1,n/2
      write(61,*) (i-1)/(n*dt),f(i)
      enddo
    close(61)


    open(61,file='out2.dat')!,status='replace',action='write')
      do i=1,n/2
      write(61,*) (i-1)/(n*dt),p(i)
      enddo
    close(61)


     !フーリエ変換後のdiffusive freqency domainのinp
    open(61,file='out3.dat')!,status='replace',action='write')
      do i=1,n/2
      write(61,*) (i-1)/(n*dt),real(c(i)),aimag(c(i))
      enddo
    close(61)

    end program f_to_d

	!!!jh-------------------------------------------------------------------
	!ファイル（データ）の長さNDを調べる
!    open(52,file=inp2,action='read')
!    nd=0
!    do
!        read(52,'(f12.0)',iostat=ios)
!        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
!         nd=nd+1
!!    enddo
 !   close(52)

 	!データの長さをnの２乗になるように決めて
    !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
!     n=2**int(log(dble(nd))/ log(2.0d0)+0.5d0)

 !    allocate(t(1:n),jh(1:n),c(1:n/2),f(1:n/2),p(1:n/2))
!    write(*,*) 'データ数nd,nの２乗数',nd,n
    
 !   if (nd>n) then
 !!   write(*,*) "nd must < n"
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



