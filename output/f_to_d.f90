!!ficticiou wave domain からdiffusive frequency domainへ

program ficticiou_to_diffusive_f
!	use const_para
		implicit none
        integer :: i,j
		integer :: n,nd,ios
		real(8) :: omega 
	!	complex(kind(0d0)), allocatable :: Gxn(:,:,:),Gyn(:,:,:),Gzn(:,:,:)
		real(8), allocatable :: jh(:)
        real(8), allocatable :: f(:),p(:),hz(:),t(:)
		complex(kind(0d0)), allocatable :: c(:)
		character(8),parameter :: inp1='hz10.d',inp2='jh.d'
		character(12),parameter :: out1='hz10as.out',out2='hz10ps.out',out3='hz10_fft.out'
        character(12),parameter :: out4='jh_as.out',out5='jh_ps.out',out6='jh_fft.out'
		integer :: plan(8)

	! FFTW3を呼び出すのに必要なヘッダーファイルを include する
    include 'fftw3.f'

    !!! Ex------------------------------------------------------------
    !ファイル（データ）の長さNDを調べる
    open(51,file=inp1,action='read')
    nd=0
    do
        read(51,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(51)

 	!データの長さをnの２乗になるように決めて
    !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
    n=2**int(log(dble(nd))/ log(2.0d0)+0.5d0)
    allocate(t(1:n),hz(1:n),c(0:n/2),f(0:n/2),p(0:n/2))
	
    open(51,file=inp1,action='read')
    do i=1,nd
    read(51,*) t(i), hz(i)
    enddo
    close(51)


 	!不足分は0パッディング
    hz(nd+1:n)=0.0d0

    !fftwの実行
    call dfftw_plan_dft_r2c_1d(plan,n,hz,c,fftw_estimate)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)


    ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
    f= abs(c)
    p= atan2(aimag(c),dble(c))

	! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
    open(61,file=out1,status='replace',action='write')
    do i=0,n/2
    write(61,*) i/(n*5.0d-4),f(i)
    enddo
    close(61)

    open(61,file=out2,status='replace',action='write')
    do i=0,n/2
    write(61,*) i/(n*5.0d-4),p(i)
    enddo
    close(61)

     !フーリエ変換後のdiffusive freqency domainのhz
    open(61,file=out3,status='replace',action='write')
    hz(1:n) = c(0:n/2)
    do i=0,n/2
    write(61,*) i/(n*5.0d-4),real(c(i)),aimag(c(i))!,aimag(c)!i, hz(i)
    enddo
    close(61)




	!!!jh-------------------------------------------------------------------
	!ファイル（データ）の長さNDを調べる
    open(52,file=inp2,action='read')
    nd=0
    do
        read(52,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(52)

 	!データの長さをnの２乗になるように決めて
    !電場exデータ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
     n=2**int(log(dble(nd))/ log(2.0d0)+0.5d0)
    
!    write(*,*) 'データ数nd,nの２乗数',nd,n
    
    if (nd>n) then
    write(*,*) "nd must < n"
    endif

    allocate(jh(1:n))
    

	open(52,file=inp2,action='read')          !!!!!!!!!!!!!!  n>ndが崩れるときがある
    do i=1,nd
    read(52,*) jh(i)
    enddo
    close(52)




 	!不足分は0パッディング
    jh(nd+1:n)=0.0d0

    !fftwの実行
    call dfftw_plan_dft_r2c_1d(plan,n,jh,c,fftw_estimate)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)



    ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
    f= abs(c)
    p= atan2(aimag(c),dble(c))

	! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
    open(61,file=out4,status='replace',action='write')
    write(61,*) f
    close(61)

    open(61,file=out5,status='replace',action='write')
    write(61,*) p
    close(61)

     !フーリエ変換後のdiffusive freqency domainのjh
        !√-2omega0/iomeagaをかける
!    jh(1:n) = sqrt(-2.0d0*omega0/(i*omega)) * c(0:n/2)
    open(61,file=out6,status='replace',action='write')
    do i=0,n/2
    write(61,*) i/(n*5.0d-4),real(c(i)),aimag(c(i))
    enddo
    close(61)

!!!グリーン関数を求める-------------------------------------------------
!	Gxn(:,:,:) = Ex / jh
!	Gyn(:,:,:) = Ey / jh
!	Gzn(:,:,:) = Ez / jh

end program ficticiou_to_diffusive_f



