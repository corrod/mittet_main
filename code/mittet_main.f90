!***************************************************
!  mittet 仮想領域法 2014.06.01
!
!  
!
!   omega0=2πf0
!   f0=1,0Hz
!   sigmawa=3.2S/m
!   x=(i-1)*dx 
!   ∂n=Σαの書き方
!   代入必要?↓↓↓
!   E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!***************************************************




!!!初期値/モデル設定,変数定義******************************************************
module const_para
    implicit none

    integer :: i,j,k
    integer, parameter :: nstep = 1000 !総タイムステップ数
    integer, parameter :: nx    = 100, ny = 100, nz = 100 !グリッド数
    integer, parameter :: x0    = 50, y0 = 50, z0 = 50  !送信源位置
    real(8), parameter :: pai   = 3.14159265358979d0 !πの値
    real(8), parameter :: fmax  = 25d0 !最大周波数
    real(8), parameter :: tau0  = 0.02!1.6d-4 !送信源出力時間
    real(8), parameter :: omega0 = 2.0d0*pai!2πf0,f0=1 !ω0
    real(8), parameter :: dt    = 5.0d-4 !タイムステップ長 s
    real(8), parameter :: dx    = 20d0,dy = 20d0,dz = 20d0!dx=1.0d-2,dy=1.0d-2,dz=1.0d-2 
    real(8)            :: sigmaxx(nx,ny,nz) !diagonal sigma x
    real(8)            :: sigmayy(nx,ny,nz) !diagonal sigma y
    real(8)            :: sigmazz(nx,ny,nz) !diagonal sigma z
    real(8), parameter :: sigmaair = 0     !空気の導電率 S/m
    real(8), parameter :: sigmafe  = 1.03d7 !鉄の導電率 S/m
    real(8), parameter :: sigmawa  = 3.2d0  !海水の導電率 S/m
    real(8), parameter :: myu0     = 1.2566370614d-6 !真空の透磁率 H/m
    real(8), parameter :: myurair  = 1.0d0        !空気の比透磁率
    real(8), parameter :: myurfe   = 4.0d3         !鉄の比透磁率
    real(8), parameter :: myurwa   = 0.999991d0     !海水の比透磁率
    real(8), parameter :: myuair   = myurair * myu0 !空気の透磁率 H/m
    real(8), parameter :: myufe    = myurfe * myu0   !鉄の透磁率 H/m
    real(8), parameter :: myuwa    = myurwa * myu0   !鉄の透磁率 H/m
    real(8), parameter :: epsi0    = 8.854d-12   !真空の誘電率F/m
    real(8), parameter :: epsirair = 1.0006d0 !空気の比誘電率
    real(8), parameter :: epsirfe  = 1.0d0     !鉄の比誘電率
    real(8), parameter :: epsirwa  = 81d0      !海水の比誘電率
    real(8), parameter :: epsiair  = epsirair * epsi0 !空気の比誘電率
    real(8), parameter :: epsife   = epsirfe * epsi0   !鉄の比誘電率
    real(8), parameter :: epsiwa   = epsirwa * epsi0   !海水の比誘電率
            end module const_para














!!!メインプログラム*************************************************************
!********************************************************************************
program main
    use const_para
    implicit none 

    integer :: istep !タイムステップ
    real(8) :: t !経過時間
    real(8) :: Jn(nstep) !gaussian
    real(8) :: Je(nstep) !電流源
    real(8) :: Jh(nstep) !磁流源
    real(8) :: sigma(nx,ny,nz),myu(nx,ny,nz)
    complex(kind(0d0)) :: Ex(nx,ny,nz)
    complex(kind(0d0)) :: Ey(nx,ny,nz)
    complex(kind(0d0)) :: Ez(nx,ny,nz)
    complex(kind(0d0)) :: Hx(nx,ny,nz)
    complex(kind(0d0)) :: Hy(nx,ny,nz)
    complex(kind(0d0)) :: Hz(nx,ny,nz)
    open(13,file='hz(x0,y0,z0+10).d') !iran
    open(14,file='hz(x0,y0,z0+20).d') !iran
    open(15,file='Jh(istep).d')!iran
    
    t=0d0!開始時間---------------------------------------

!   Ex,Ey,Ez,Hx,Hy,Hzの初期化
    do k=1,nz
        do j=1,ny
            do i=1,nx
                Ex(i,j,k)=0d0
                Ey(i,j,k)=0d0
                Ez(i,j,k)=0d0
                Hx(i,j,k)=0d0
                Hy(i,j,k)=0d0
                Hz(i,j,k)=0d0
            enddo 
        enddo
    enddo

    !モデルの読み込み
    call model(sigma,myu)

    !cmax,cminの計算 dt,dx,dy,dzの設定 
    call set_d_txyz

    t = 0d0  !時間の初期化

    do istep = 1, nstep !*反復計算開始------------------------

    !入力波源の設定
!   call gaussianpulse(istep,t,Ie,Mh)  
    call gaussian(istep,t,Je,Jh,sigma,myu)
    !電場計算
    call EXFIELD(istep,t,Je,Ex,Hy,Hz,sigma)
    call EYFIELD(istep,t,Je,Ey,Hz,Hx,sigma)
    call EZFIELD(istep,t,Je,Ez,Hx,Hy,sigma)

     !境界条件 E
!   call E_PML(Ex,Ey,Ez,HX,Hy,Hz,sigma)

    t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算
    call HXFIELD(istep,t,Jh,Hx,Ey,Ez,myu)
    call HYFIELD(istep,t,Jh,Hy,Ex,Ez,myu)
    call HZFIELD(istep,t,Jh,Hz,Ex,Ey,myu)

    !境界条件 H
!   call H_BoundaryCondition(Ex,Ey,Ez,Hx,Hy,Hz,sigma,myu)

    t = t + dt*0.5d0 !時間の更新---------------------------

    !アウトプットE-field、H-field
    call output_EH(istep,t,Ex,Ey,Ez,Hx,Hy,Hz)

    enddo !*反復計算終了


    !グリーン関数の導出
    !call green()

    !高速フーリエ変換
    !call fft()

    close(13)
    close(14)
    close(15)
            endprogram main





















!!!dt,dx,dy,dzの設定cmax,cminの計算************************************************
subroutine set_d_txyz
    use const_para
    implicit none

	real(8) :: cwa, cfe, cair
	real(8) :: cmax
    real(8) :: t_max
	real(8) :: dt_max

    !媒質中の伝播速度計算
    cwa = sqrt(2.0d0*omega0/myuwa/sigmawa)
    cfe = sqrt(2.0d0*omega0/myufe/sigmafe)
    cair = 1.0d0 / sqrt(myuair*epsiair)
    write(*,'(a,2e12.4)') 'cwa,cfe', cwa,cfe !海水の伝播速度出力
    write(*,'(a,e12.4)') '伝播距離c*t', cwa*dt*nstep

  !  cmax = max(cwa,cfe,cair)!   最大伝播速度cmax計算
  !  write(*,*) 'cmax=',cmax

    t_max = dt*nstep !   計測時間t_max
    write(*,'(a,e12.4)') '計測時間t_max',t_max
    dt_max = (2.0d0*dx) / ((3**0.5)*3.14159d0*cwa) !  
    write(*,'(a,e12.4)') 'dt_max', dt_max !最大タイムステップ長の出力
    write(*,'(a,e12.4)') '１辺の長さm',dx*nx 
    write(*,'(a,e12.4)') 'dt',dt !タイムステップの出力
    write(*,'(a,3e12.4)') 'dx,dy,dz',dx,dy,dz !グリッドサイズの出力
    write(*,'(a,e12.4)') 'tau0',tau0 !ソースの継続時間
    write(*,'(a,3i5)') 'nx,ny,nz',nx,ny,nz !グリッド数
    write(*,'(a,i5)') 'nstep',nstep
    write(*,'(a,e12.4)') 'omega0',omega0
    write(*,*) 'sigmawa,myuwa',sigmawa,myuwa
            endsubroutine set_d_txyz






















!!!送信源the first derivative of Gaussian********************************
subroutine gaussian(istep,t,Je,Jh,sigma,myu)
    use const_para
    implicit none
    
    integer, intent(in)  :: istep
    real(8), intent(in)  :: t
    real(8)              :: Jn(nstep)!gaussian
    real(8), intent(out) :: Je(nstep)
    real(8), intent(out) :: Jh(nstep)
    real(8), intent(in)  :: sigma(nx,ny,nz)
    real(8), intent(in)  :: myu(nx,ny,nz)
    real(8)              :: etaxx(x0,y0,z0)
    real(8), parameter   :: t0 = pai/fmax
    real(8), parameter   :: beta = pai*(fmax**2)  

    Jn(istep) = -(2.0d0*beta*(t-t0)*sqrt(beta/pai))*exp(-beta*(t-t0)**2)

    !電場ソースの設定
    etaxx(x0,y0,z0) = (2.0d0*omega0)/sigma(x0,y0,z0)
    Je(istep) = dt*etaxx(x0,y0,z0)*Jn(istep)/dx/dy/dz

    !磁場ソースの設定
    Jh(istep) = Jn(istep)*dt/myu(x0,y0,z0)/dx/dy/dz

    write(12,*) 't,Jn,Je,Jh',t,Jn(istep) 
    write(15,*) t,Jh(istep)
        endsubroutine gaussian





















!!!モデル設定********************************************************************
subroutine model(sigma,myu)
    use const_para
    implicit none

	real(8),intent(out) :: sigma(nx,ny,nz)
	real(8),intent(out) :: myu(nx,ny,nz)

    !海水一様モデル
    sigma(1:nx,1:ny,1:nz) = sigmawa
    myu(1:nx,1:ny,1:nz) = myuwa

    !海水
    !sigma()=sigmawa
    !myu()=myuwa

    !鉄板
    !sigma()=sigmafe
    !myu()=myufe
            endsubroutine model





















!!!アウトプットE-field,H-field******************************************************
subroutine output_EH(istep,t,Ex,Ey,Ez,Hx,Hy,Hz)
    use const_para
    implicit none

    integer :: l
    integer,intent(in) :: istep
    real(8), intent(in) :: t
    complex(kind(0d0)), intent(in) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
    character(5) :: name
 !   open(1,file='Exfield.d', position='append')
 !   open(2,file='Eyfield.d', position='append')
 !   open(3,file='Ezfield.d', position='append')
  !  open(4,file='Hxfield.d', position='append')
  !  open(5,file='Hyfield.d', position='append')
  !  open(6,file='Hzfield.d', position='append')

!    do k=1,nz
!         do j=1,ny
!             do i=1,nx
!               write(1,'(i5,3i4,e12.4)') istep,i,j,k,real(Ex(i,j,k))
!               write(2,'(i5,3i4,e12.4)') istep,i,j,k,real(Ey(i,j,k))
!               write(3,'(i5,3i4,e12.4)') istep,i,j,k,real(Ez(i,j,k))
!               write(4,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!               write(5,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!               write(6,'(i5,3i4,e12.4)') istep,i,j,k,real(Hx(i,j,k))
!            enddo
!          enddo
!    enddo

 !   close(1)
  !  close(2)
   ! close(3)
!    close(4)
!    close(5)
!    close(6)

    write(13,*) t, real(hz(x0,y0,z0+10)) !iran
    write(14,*) t, real(hz(x0,y0,z0+20)) !iran

!--------シェル用出力------
    if (mod(istep,20)==0) then
   l=10000+istep/20
    write(name,"(I5)") l
    open(7,file=name//".d")
        do j=1,ny
            do i=1,nx
                 write(7,*) t,i,j,real(Ex(i,j,z0))
            enddo
        enddo
    close(7)
    endif
            endsubroutine output_EH














!!!仮想領域での電磁場の計算***********************************************************

!電場計算*****************************************************************************
!Ex-field
subroutine EXFIELD(istep,t,Je,Ex,Hy,Hz,sigma)
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)), intent(in) :: sigmaxx(nx,ny,nz)
    real(8), intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    real(8),intent(in)  :: sigma(nx,ny,nz)
    real(8) :: etaxx(nx,ny,nz)
    real(8) :: CEXLY(nx,ny,nz)
    real(8) :: CEXLZ(nx,ny,nz)

    !係数の設定    
    do k=1,nz
         do j=1,ny
              do i=1,nx
                    etaxx(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmaxx(i,j,k)
                    CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
                    CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
              enddo
         enddo
    enddo

    !電場ソースの設定
    !Je(istep) = dt * etaxx(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Ex(i,j,k) = Ex(i,j,k)&
                            + CEXLY(i,j,k) * (Hz(i,j,k) - Hz(i,j-1,k))&
                            + CEXLZ(i,j,k) * (Hy(i,j,k) - Hy(i,j,k-1))
            enddo
        enddo
    enddo

    !ソース項
    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)

   !ソース位置、ソースから離れた位置でのExの磁場の時間分布
   write(9,*) 't,ex(x0+10,y0,z0)',t,real(ex(x0+10,y0,z0))
            endsubroutine EXFIELD


!Ey-field-------------------------------------------------------------
subroutine EYFIELD(istep,t,Je,Ey,Hz,Hx,sigma)
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)), intent(in) :: sigmayy(nx,ny,nz)
    real(8), intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8) :: etayy(nx,ny,nz)
    real(8) :: CEYLZ(nx,ny,nz)
    real(8) :: CEYLX(nx,ny,nz)

    !係数の設定    
    do k=1,nz
        do j=1,ny
           do i=1,nx
              etayy(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmayy(i,j,k)
              CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
              CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
           enddo
        enddo
    enddo

    !電場ソースの設定
    !  Je(istep) = dt * etayy(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                Ey(i,j,k) = Ey(i,j,k)&
                            + CEYLZ(i,j,k) * (Hx(i,j,k) - Hx(i,j,k-1))&
                            + CEYLX(i,j,k) * (Hz(i,j,k) - Hz(i-1,j,k))
            enddo
        enddo
    enddo

    !ソース項
    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)
            endsubroutine EYFIELD



!Ez-field-------------------------------------------------------------------
subroutine EZFIELD(istep,t,Je,Ez,Hx,Hy,sigma)
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)),intent(in) :: sigmazz(nx,ny,nz)
    real(8),intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    real(8),intent(in)  :: sigma(nx,ny,nz)
    real(8) :: etazz(nx,ny,nz)
    real(8) :: CEZLX(nx,ny,nz)
    real(8) :: CEZLY(nx,ny,nz)

    !係数の設定    
    do k=1,nz
        do j=1,ny
            do i=1,nx
                etazz(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmazz(i,j,k)
                CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
                CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
            enddo
        enddo
    enddo

    !電場ソースの設定
 !   Je(istep) = dt * etazz(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                Ez(i,j,k) = Ez(i,j,k)&
                            + CEZLX(i,j,k) * (Hy(i,j,k) - Hy(i-1,j,k))&
                            + CEZLY(i,j,k) * (Hx(i,j,k) - Hx(i,j-1,k))
            enddo
        enddo
    enddo

    !ソース項
!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)
            endsubroutine EZFIELD













!磁場計算*********************************************************************
!Hx-field
subroutine  HXFIELD(istep,t,Jh,Hx,Ey,Ez,myu)!myu追加
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz) 
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8) :: CHXLY(nx,ny,nz)
    real(8) :: CHXLZ(nx,ny,nz)

    !係数の設定    
    do k=1,nz
       do j=1,ny
           do i=1,nx
              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
          enddo
      enddo
    enddo

    !磁場ソースの設定
   ! Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz
   

    !波動伝播計
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Hx(i,j,k) = Hx(i,j,k)&
                        + CHXLY(i,j,k) * (Ez(i,j+1,k) - Ez(i,j,k))&
                        + CHXLZ(i,j,k) * (Ey(i,j,k+1) - Ey(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)
            endsubroutine HXFIELD




!Hy-field-------------------------------------------------------------
subroutine HYFIELD(istep,t,Jh,Hy,Ex,Ez,myu) !myu追加
    use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz)
    real(8) :: CHYLZ(nx,ny,nz)
    real(8) :: CHYLX(nx,ny,nz)

    !係数の設定    
    do k=1,nz
        do j=1,ny
            do i=1,nx
             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
          enddo
      enddo
    enddo

    !磁場ソースの設定
  !  Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Hy(i,j,k) = Hy(i,j,k)&
                        + CHYLZ(i,j,k) * (Ex(i,j,k+1) - Ex(i,j,k))&
                        + CHYLX(i,j,k) * (Ez(i+1,j,k) - Ez(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
!   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)
            endsubroutine HYFIELD




!Hz-field--------------------------------------------------------------------------
subroutine HZFIELD(istep,t,Jh,Hz,Ex,Ey,myu) !myu追加
    use const_para
    implicit none

    integer :: l !スクリプト用
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    real(8) :: CHZLX(nx,ny,nz)
    real(8) :: CHZLY(nx,ny,nz)
    character(8) :: name

    !係数の設定    
    do k=1,nz
       do j=1,ny
         do i=1,nx
            CHZLX(i,j,k) = - dt / myu(i,j,k) / dz
            CHZLY(i,j,k) = dt / myu(i,j,k) / dx
        enddo
      enddo
    enddo

    !磁場ソースの設定
   ! Jh(istep) = Jn(istep) * dt /myu(x0,y0,z0) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                 Hz(i,j,k) = Hz(i,j,k)&
                        + CHZLX(i,j,k) * (Ey(i+1,j,k) - Ey(i,j,k))&
                        + CHZLY(i,j,k) * (Ex(i,j+1,k) - Ex(i,j,k))
            enddo
        enddo
    enddo

    !ソース項
    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)


    !ソース位置、ソースから離れた位置でのhzの磁場の時間分布
   write(8,*) 't,jh,hz(x0,y0,z0+10)',t,Jh(istep),real(hz(x0,y0,z0+10))

!--------シェル用出力------
!   if (mod(istep,5)==0) then
!    l=10000+istep/5
!    write(name,"(I5)") l
!    open(9,file=name//".d")
!    do j=1,ny
!        do i=1,nx
!            write(2,*) i,j,real(hz(i,j,z0))
!        enddo
!    enddo
!    close(9)
!   endif
            endsubroutine HZFIELD



























!Convolutional PML_E **********************************************************************************
subroutine cpml
    use const_para
    implicit none
 
    integer,parameter :: m=4,ma=4
    integer,parameter :: nxpml1=5,nypml1=5,nzpml1=5
    real(8) :: sigma_max !!! 導出法確認
    real(8) :: kappa_max !!! 整数？
    real(8) :: a_max     !!!
    real(8) :: sigma(nx,ny,nz)
    real(8) :: sigma_x(nx),sigma_y(ny),sigma_z(nz)
    real(8) :: ca_x(nx,ny,nz),ca_y(nx,ny,nz),ca_z(nx,ny,nz)
    real(8) :: cb_x(nx,ny,nz),cb_y(nx,ny,nz),cb_z(nx,ny,nz)
    real(8) :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8) :: kedx(nx),kedy(ny),kedz(nz)
    real(8) :: epsi(nx,ny,nz)
    real(8) :: a_x(nx),a_y(ny),a_z(nz)
    real(8) :: be_x(nx),be_y(ny),be_z(nz)
    real(8) :: ce_x(nx),ce_y(ny),ce_z(nz)
    complex(kind(0d0)) :: psi_ezx1(nx,ny,nz),psi_eyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_exy1(nx,ny,nz),psi_ezy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_eyz1(nx,ny,nz),psi_exz1(nx,ny,nz)
    complex(kind(0d0)) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz) 
    complex(kind(0d0)) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

!係数の設定
    ca_x(i,j,k) = (1-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    ca_y(i,j,k) = (1-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    ca_z(i,j,k) = (1-sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))

    cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
    cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1+sigma(i,j,k)*dt/(2.0d0*epsi(i,j,k)))

    kedx(i) = kappa_x(i)*dx  !kedy=kappa(jdy)dy
    kedy(j) = kappa_y(j)*dy  !!!
    kedz(k) = kappa_z(k)*dz  !!!
   
    sigma_x(i) = sigma_max*((nxpml1-i)/(nxpml1-1))**m
    sigma_y(j) = sigma_max*((nypml1-j)/(nypml1-1))**m
    sigma_z(k) = sigma_max*((nzpml1-k)/(nzpml1-1))**m

    kappa_x(i) = 1 + (kappa_max-1)*((nxpml1-j)/(nxpml1-1))**m
    kappa_y(j) = 1 + (kappa_max-1)*((nypml1-j)/(nypml1-1))**m
    kappa_z(k) = 1 + (kappa_max-1)*((nzpml1-j)/(nzpml1-1))**m
    
    a_x(j) = a_max*((i-1)/(nxpml1-1))**ma
    a_y(j) = a_max*((j-1)/(nypml1-1))**ma
    a_z(j) = a_max*((k-1)/(nzpml1-1))**ma

    be_x(i) = exp((sigma_x(i)/kappa_x(i)+a_x(i))*dt/epsi0)
    be_y(j) = exp((sigma_y(j)/kappa_y(j)+a_y(j))*dt/epsi0)
    be_z(k) = exp((sigma_z(k)/kappa_z(k)+a_z(k))*dt/epsi0)

    ce_x(i) = sigma_x(i)*(be_x(i)-1) / (sigma_x(i)+kappa_x(i)*a_x(i)) / kappa_x(i)
    ce_y(j) = sigma_y(j)*(be_y(j)-1) / (sigma_y(j)+kappa_y(j)*a_y(j)) / kappa_y(j)
    ce_z(k) = sigma_z(k)*(be_z(k)-1) / (sigma_z(k)+kappa_z(k)*a_z(k)) / kappa_z(k)

!field-update loop
    !x-update
    do  k=2,nz-1
        do  j=2,ny-1
            do  i=1,nx-1
                ex(i,j,k) = ca_x(i,j,k)*ex(i,j,k)&
                            +cb_x(i,j,k)*((hz(i,j,k)-hz(i,j-1,k)) / kedy(j) &
                                            - (hy(i,j,k)-hy(i,j,k-1))/ kedz(k))
            enddo
        enddo
    enddo                    
    !y-update
    do  k=2,nz-1
        do  j=1,ny-1
            do i=2,nx-1
                ey(i,j,k) = ca_y(i,j,k)*ey(i,j,k)&
                            +cb_y(i,j,k)*((hx(i,j,k)-hx(i,j-1,k)) / kedy(j) &  !インデント
                                           - (hz(i,j,k)-hz(i,j,k-1))/ kedz(k)) !kedyとか
            enddo
        enddo
    enddo
    !z-update
     do  k=1,nz-1
        do  j=2,ny-1
            do  i=2,nx-1
                ez(i,j,k) = ca_z(i,j,k)*ez(i,j,k)&
                            +cb_z(i,j,k)*((hy(i,j,k)-hy(i,j-1,k)) / kedy(j) &  !インデント
                                           - (hx(i,j,k)-hx(i,j,k-1))/ kedz(k)) !kedyとか
            enddo
        enddo
     enddo

!psi-update
    !x-PML loop
    do k =1,nz-1
        do j=1,ny-1
            do i=2,nxpml1
                psi_ezx1(i,j,k) = be_x(j)*psi_ezx1(i,j,k)&
                                +ce_x(j)*(hy(i,j,k)-hy(i-1,j,k)) / dx
                psi_eyx1(i,j,k) = be_x(j)*psi_eyx1(i,j,k)&
                                +ce_x(j)*(hz(i,j,k)-hz(i-1,j,k)) / dx
                ez(i,j,k) = ez(i,j,k)+cb_z(i,j,k)*psi_ezx1(i,j,k)
                ey(i,j,k) = ey(i,j,k)-cb_y(i,j,k)*psi_eyx1(i,j,k)
            enddo
        enddo
    enddo

    !y-PML loop
    do k =1,nz-1
        do j=2,nypml1
            do i=1,nx-1
                psi_exy1(i,j,k) = be_y(j)*psi_exy1(i,j,k)&
                                +ce_y(j)*(hz(i,j,k)-hz(i,j-1,k)) / dy
                psi_ezy1(i,j,k) = be_y(j)*psi_ezy1(i,j,k)&
                                +ce_y(j)*(hx(i,j,k)-hx(i,j-1,k)) / dy
                ex(i,j,k) = ex(i,j,k)+cb_x(i,j,k)*psi_exy1(i,j,k)
                ez(i,j,k) = ez(i,j,k)-cb_z(i,j,k)*psi_ezy1(i,j,k)
            enddo
        enddo
    enddo

    !z-PML loop
    do k =2,nzpml1
        do j=1,ny-1
            do i=1,nx-1
                psi_eyz1(i,j,k) = be_z(j)*psi_eyz1(i,j,k)&
                                +ce_z(j)*(hx(i,j,k)-hx(i,j,k-1)) / dz
                psi_exz1(i,j,k) = be_z(j)*psi_exz1(i,j,k)&
                                +ce_z(j)*(hy(i,j,k)-hy(i,j,k-1)) / dz
                ey(i,j,k) = ey(i,j,k)+cb_y(i,j,k)*psi_eyz1(i,j,k)
                ex(i,j,k) = ex(i,j,k)-cb_x(i,j,k)*psi_exz1(i,j,k)
            enddo
        enddo
    enddo
        endsubroutine cpml








!Convolutional PML_H ********************************************************************************************
subroutine cpml_H
    use const_para
    implicit none

    integer,parameter :: m=4, ma=4
    integer :: nxpml1=5, nypml1=5, nzpml1=5 !pmlの厚さ
    real(8),parameter :: lnR0= -100d0  !ln|R(0)|
    real(8) :: sigma_max !!!
    real(8) :: kappa_max !!!
    real(8) :: a_max     !!!
    real(8) :: myu(nx,ny,nz),epsi(nx,ny,nz)
    real(8) :: da_x(nx,ny,nz),da_y(nx,ny,nz),da_z(nx,ny,nz)
    real(8) :: db_x(nx,ny,nz),db_y(nx,ny,nz),db_z(nx,ny,nz)
    real(8) :: khdx(nx),khdy(ny),khdz(nz)
    real(8) :: bh_x(nx),bh_y(ny),bh_z(dz)
    real(8) :: ch_x(nx),ch_y(ny),ch_z(nz)
    real(8) :: sigma(nx,ny,nz)
    real(8) :: sigma_x(nx),sigma_y(ny),sigma_z(nz)
    real(8) :: kappa_x(nx),kappa_y(ny),kappa_z(nz)
    real(8) :: a_x(nx),a_y(ny),a_z(nz)
    complex(kind(0d0)) :: psi_hzx1(nx,ny,nz),psi_hyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_hxy1(nx,ny,nz),psi_hzy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_hyz1(nx,ny,nz),psi_hxz1(nx,ny,nz)
    complex(kind(0d0)) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)
    complex(kind(0d0)) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)

!係数の設定
    da_x(i,j,k) = (1-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sigma=σ*
    da_y(i,j,k) = (1-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
    da_z(i,j,k) = (1-(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))

    db_x(i,j,k) = (dt/myu(i,j,k)) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_y(i,j,k) = (dt/myu(i,j,k)) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    db_z(i,j,k) = (dt/myu(i,j,k)) / (1+(sigma(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
      
!!!    sigma_max = -(m+1)*lnR0 / (2.0d0*(sqrt(myu/epsi))*nxpml1*dx)  !ln(R(0));反射係数!!!

    sigma_x(i) = sigma_max* ((nxpml1-i-1/2)/(nxpml1-1)) **m
    sigma_y(j) = sigma_max* ((nypml1-j-1/2)/(nypml1-1)) **m
    sigma_z(k) = sigma_max* ((nzpml1-j-1/2)/(nzpml1-1)) **m

!!!    kappa_max =    !!!導出要確認
    kappa_x(i) = 1 + (kappa_max-1)*((nxpml1-i-1/2)/(nxpml1-1)) **m
    kappa_y(j) = 1 + (kappa_max-1)*((nypml1-j-1/2)/(nypml1-1)) **m
    kappa_z(k) = 1 + (kappa_max-1)*((nzpml1-k-1/2)/(nzpml1-1)) **m

!!!    a_max =         !!導出要確認
    a_x(i) = a_max* ((i-1/2)/(nxpml1-1))**ma
    a_y(j) = a_max* ((j-1/2)/(nypml1-1))**ma
    a_z(k) = a_max* ((k-1/2)/(nzpml1-1))**ma

    khdx(i) = kappa_x((i-1/2)*dx)*dx
    khdy(i) = kappa_y((j-1/2)*dy)*dy
    khdz(i) = kappa_z((k-1/2)*dy)*dz

    bh_x(i) = exp((sigma_x(i)/kappa_x(i)+a_x(i)) *dt/epsi0)  !!カッコ位置確認
    bh_y(j) = exp((sigma_y(j)/kappa_y(j)+a_y(j)) *dt/epsi0)
    bh_z(k) = exp((sigma_z(k)/kappa_z(k)+a_z(k)) *dt/epsi0)

    ch_x(i) = sigma_x(j)*(bh_x(i)-1) / (sigma_x(i) + kappa_x(i)*a_x(i)) / kappa_x(i)
    ch_y(j) = sigma_y(j)*(bh_y(j)-1) / (sigma_y(j) + kappa_y(j)*a_y(j)) / kappa_y(j)
    ch_z(j) = sigma_z(k)*(bh_z(k)-1) / (sigma_z(k) + kappa_z(k)*a_z(k)) / kappa_z(k)

!field update loop 
 !Hx
 do k=1,nz-1
        do j=1,ny-1
            do i=2,nx-1
                hx(i,j,k) = da_x(i,j,k)*hx(i,j,k)&
                            -db_x(i,j,k)*((ez(i,j+1,k)-ez(i,j,k)) /khdy(j)-((ey(i,j,k+1)-ey(i,j,k))/khdz(k)))
                            enddo
                                enddo
                                    enddo
 !Hy                                   
 do k=1,nz-1
        do j=2,ny-1
            do i=1,nx-1
                hy(i,j,k) = da_y(i,j,k)*hy(i,j,k)&
                            -db_y(i,j,k)*((ex(i,j,k+1)-ex(i,j,k)) /khdz(k)-((ez(i+1,j,k)-ez(i,j,k))/khdx(i)))
                            enddo
                                enddo
                                    enddo
 !Hz                                    
 do k=2,nz-1
        do j=1,ny-1
            do i=1,nx-1
                hz(i,j,k) = da_z(i,j,k)*hz(i,j,k)&
                            -db_z(i,j,k)*((ey(i+1,j,k)-ey(i,j,k)) /khdz(i)-((ex(i,j+1,k)-ex(i,j,k))/khdx(j)))
                            enddo
                                enddo
                                    enddo

!psi update
!x-PML loop
    do k=1,nz-1
        do j=1,ny-1
            do i=1,nxpml1-1
                psi_hzx1(i,j,k) = bh_x(i)*psi_hzx1(i,j,k)&
                                   +ch_x(i)*(ey(i+1,j,k)-ey(i,j,k)) / dx
                psi_hyx1(i,j,k) = bh_z(j)*psi_hxz1(i,j,k)&
                                   +ch_x(k)*(ez(i+1,j,k)-ez(i,j,k)) / dx
                hz(i,j,k) = hz(i,j,k)-db_z(i,j,k)*psi_hzx1(i,j,k)
                hy(i,j,k) = hy(i,j,k)+db_y(i,j,k)*psi_hyx1(i,j,k)
                   enddo
                        enddo
                            enddo

!y-PML loop
    do k=1,nz-1
        do j=1,nypml1-1
            do i=1,nx-1
                psi_hxy1(i,j,k) = bh_y(j)*psi_hxy1(i,j,k)&
                                +ch_y(j)*(ez(i,j+1,k)-ez(i,j,k)) / dy
                psi_hzy1(i,j,k) = bh_y(j)*psi_hzy1(i,j,k)&
                                +ch_y(j)*(ex(i,j+1,k)-ex(i,j,k)) / dy
                hx(i,j,k) = hx(i,j,k)-db_x(i,j,k)*psi_hxy1(i,j,k)
                hz(i,j,k) = hz(i,j,k)+db_z(i,j,k)*psi_hzy1(i,j,k)
                  enddo
                       enddo
                           enddo
!z-PML loop
    do k=1,nzpml1-1
        do j=1,ny-1
            do i=1,nx-1
                psi_hyz1(i,j,k) = bh_z(k)*psi_hyz1(i,j,k)&
                                +ch_z(k)*(ex(i,j,k+1)-ex(i,j,k)) / dz
                psi_hxz1(i,j,k) = bh_y(j)*psi_hzy1(i,j,k)&
                                +ch_z(k)*(ey(i,j,k+1)-ey(i,j,k)) / dz
                hy(i,j,k) = hy(i,j,k)-db_y(i,j,k)*psi_hyz1(i,j,k)
                hx(i,j,k) = hx(i,j,k)+db_x(i,j,k)*psi_hxz1(i,j,k)
                   enddo
                        enddo
                            enddo
            endsubroutine cpml_H

































!!!境界条件 E_PML***************************************************************************
!係数の導出がまだ
!cとは？
!*******************************************************************************************
subroutine E_PML(Ex,Ey,Ez,Hx,Hy,Hz,c)
    use const_para
    implicit none

    integer :: I0,I1,J0,J1,K0,K1
    integer :: L1,L2,L3,L
    integer :: LPMLII(6,2)
    integer :: LPMLJJ(6,2)
    integer :: LPMLKK(6,2)
    integer :: LPMLST(6)
    integer,parameter :: L_max=6 !層数
    integer :: M=4 !導電率の分布を与える次数
    real(8) :: sigma_max
    real(8) :: sigmax(L_max),sigmay(L_max),sigmaz(L_max)
    real(8) :: epsi !sigma/2omega0で置き換えれるかも
    real(8) :: c
    real(8) :: R !反射係数
    complex(kind(0d0)),intent(inout) :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
    complex(kind(0d0)),intent(in)    :: Hx(nx,ny,nz), Hy(nx,ny,nz), Hz(nx,ny,nz)
    complex(kind(0d0)) :: Exy(nz),Exz(ny),Eyx(nz),Eyz(nx),Ezx(ny),Ezy(nx)
    complex(kind(0d0)) :: CEXY(ny),CEXYL(ny),CEXZ(nz),CEXZL(nz)!**配列の大きさ適当
    complex(kind(0d0)) :: CEYX(nx),CEYXL(nx),CEYZ(nz),CEYZL(nz)!**配列の大きさ適当
    complex(kind(0d0)) :: CEZX(nx),CEZXL(nx),CEZY(ny),CEZYL(ny)!**配列の大きさ適当

    sigma_max = - (((M+1)*epsi*c) / (2.0d0*L_max*dx)) * R   
     do i=1,L_max
          sigmax(i) = sigma_max * ((i-0.5d0)/L_max)
     enddo
     do j=1,L_max
          sigmay(j) = sigma_max * ((j-0.5d0)/L_max)
     enddo
     do k=1,L_max
          sigmaz(k) = sigma_max * ((k-0.5d0)/L_max)
     enddo


   !mittet sigma=2*omega0*epsi PML
    CEXY(j)  = (1-omega0*dt) / (1+omega0*dt)
    CEXYL(j) = 2.0d0*omega0*dt / (sigmay(j)*(1+dt*omega0))
    CEXZ(k)  = (1-omega0*dt) / (1+omega0*dt)
    CEXZL(k) = 2.0d0*omega0*dt / (sigmaz(k)*(1+dt*omega0))

    CEYX(i)  = (1-omega0*dt) / (1+omega0*dt)
    CEYXL(i) = 2.0d0*omega0*dt / (sigmax(i)*(1+dt*omega0))
    CEYZ(k)  = (1-omega0*dt) / (1+omega0*dt)
    CEYZL(k) = 2.0d0*omega0*dt / (sigmaz(k)*(1+dt*omega0))

    CEZX(i)  = (1-omega0*dt) / (1+omega0*dt)
    CEZXL(i) = 2.0d0*omega0*dt / (sigmax(i)*(1+dt*omega0))
    CEZY(j)  = (1-omega0*dt) / (1+omega0*dt)
    CEZYL(j) = 2.0d0*omega0*dt / (sigmay(j)*(1+dt*omega0))

    !normal PML 
    !CEXY(j) = (1-sigma(j)*dt/(2.0d0*epsi)) / (1+sigma(j)*dt/(2.0d0*epsi))
    !CEXYL(j) = (dt/epsi) / (1+sigma(j)*dt/(2.0d0*epsi)) / dy
    !CEXZ(k) =(1-sigma(j)*dt/(2.0d0*epsi)) / (1+sigma(j)*dt/(2.0d0*epsi))
    !CEXZL(k)=(dt/epsi) / (1+sigma(k)*dt/(2.0d0*epsi)) / dz

    !CEYX(i) =(1-sigma(i)*dt/(2.0d0*epsi)) / (1+sigma(i)*dt/(2.0d0*epsi))
    !CEYXL(i) =(dt/epsi) / (1+sigma(i)*dt/(2.0d0*epsi)) / dx
    !CEYZ(k) =(1-sigma(k)*dt/(2.0d0*epsi)) / (1+sigma(k)*dt/(2.0d0*epsi))
    !CEYZL(k) =(dt/epsi) / (1+sigma(k)*dt/(2.0d0*epsi)) / dz

    !CEZX(i) =(1-sigma(i)*dt/(2.0d0*epsi)) / (1+sigma(i)*dt/(2.0d0*epsi))
    !CEZXL(i) =(dt/epsi) / (1+sigma(i)*dt/(2.0d0*epsi)) / dx
    !CEZY(j) =(1-sigma(j)*dt/(2.0d0*epsi)) / (1+sigma(j)*dt/(2.0d0*epsi))
    !CEZYL(j) =(dt/epsi) / (1+sigma(j)*dt/(2.0d0*epsi)) / dy


     do L=1,6
      I0 = LPMLII(L,1)
      I1 = LPMLII(L,2)
      J0 = LPMLJJ(L,1)
      J1 = LPMLJJ(L,2)
      K0 = LPMLKK(L,1)
      K1 = LPMLKK(L,2)

    L1=LPMLST(L)
    do i=I0,I1-1
        do j=J0+1,J1-1
            do k=K0+1,K1-1
                Exy(L1) = CEXY(J)*Exy(L1)&
                         +CEXYL(j)*(Hz(i,j,k)-Hz(i,j-1,k))
                Exz(L1) = CEXZ(k)*Exz(L1)&
                         +CEXZL(k)*(Hy(i,j,k-1)-Hy(i,j,k))
                Ex(i,j,k) = Exy(L1)+Exz(L1)
                L1=L1+1
            enddo
        enddo
    enddo
  
    L2=LPMLST(L)
    do i=I0+1,I1-1
        do j=J0,J1-1
            do k=K0+1,K1-1
                Eyx(L2) = CEYX(I)*EYX(L2)&
                          +CEYXL(i)*(HZ(i-1,j,k)-Hz(i,j,k))
                Eyz(L2) = CEYZ(k)*Eyz(L2)&
                          +CEYZL(k)*(Hx(i,j,k)-Hx(i,j,k-1))
                Ey(i,j,k) = Eyx(L2)+Eyz(L2)
                L2=L2+1
            enddo
        enddo
    enddo

   L3=LPMLST(L)
    do i=I0+1,I1-1
        do j=J0+1,J1-1
            do k=K0+1,K1-1
                Ezx(L3) = CEZX(I)*Ezx(L3)&
                          +CEZXL(i)*(HY(i,j,k)-Hy(i-1,j,k))
                Ezy(L3) = CEZY(j)*Ezy(L3)&
                          +CEZYL(j)*(Hx(i,j-1,k)-Hx(i,j,k))
                Ez(i,j,k) = Ezx(L3)+Ezy(L3)
                L3=L3+1
            enddo
          enddo
       enddo
   
    enddo
                endsubroutine E_PML













!!!境界条件H PML*********************************************************************************
!subroutine H_BoundaryCondition(Ex,Ey,Ez,Hx,Hy,Hz,sigma,myu)
!    use const_para
!    implicit none
!    integer :: I0,I1,J0,J1,K0,K1
!    integer :: L1,L2,L3,L
!    integer :: LPMLII(6,2)
 !   integer :: LPMLJJ(6,2)
 !   integer :: LPMLKK(6,2)
 !   integer :: LPMLST(6)
 !   real(8), intent(in) :: sigma, myu
  !  complex(kind(0d0)),intent(in)    :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
  !  complex(kind(0d0)),intent(inout) :: Hx(nx,ny,nz), Hy(nx,ny,nz), Hz(nx,ny,nz)
 !   complex(kind(0d0)) :: Hxy(nz),Hxz(ny),Hyx(nz),Hyz(nx),Hzx(ny),Hzy(nx)
   ! complex(kind(0d0)) :: CHXY(ny),CHXYL(ny),CHXZ(nz),CHXZL(nz)!**配列の大きさ適当
    !complex(kind(0d0)) :: CHYX(nx),CHYXL(nx),CHYZ(nz),CHYZL(nz)!**配列の大きさ適当
    !complex(kind(0d0)) :: CHZX(nx),CHZXL(nx),CHZY(ny),CHZYL(ny)!**配列の大きさ適当
   
   !mittet sigma=2*omega0*epsi PML
    !CHXY(j)  = 
    !CHXYL(j) = 
    !CHXZ(k)  = 
    !CHXZL(k) = 

    !CHYX(i)  = 
    !CHYXL(i) = 
    !CHYZ(k)  = 
    !CHYZL(k) = 

    !CHZX(i)  = 
    !CHZXL(i) = 
    !CHZY(j)  = 
    !CHZYL(j) =

     !normal PML **インデックス要確認 
   ! CHXY(j)  = (1-sigma(j)*dt/(2.0*epsi)) / (1+sigma(j)*dt/(2.0d0epsi))
    !CHXYL(j) = -(dt/myu) / (1+sigma(j)*dt/(2.0d0*epsi))/dx 
    !CHXZ(k)  = (1-sigma(k)*dt/(2.0*epsi)) / (1+sigma(k)*dt/(2.0d0epsi))
    !CHXZL(k) = -(dt/myu) / (1+sigma(k)*dt/(2.0d0*epsi))/dx 

    !CHYX(i)  = (1-sigma(i)*dt/(2.0*epsi)) / (1+sigma(i)*dt/(2.0d0epsi))
    !CHYXL(i) = -(dt/myu) / (1+sigma(i)*dt/(2.0d0*epsi))/dx 
    !CHYZ(k)  = (1-sigma(k)*dt/(2.0*epsi)) / (1+sigma(k)*dt/(2.0d0epsi))
    !CHYZL(k) = -(dt/myu) / (1+sigma(k)*dt/(2.0d0*epsi))/dx 

    !CHZX(i)  = (1-sigma(i)*dt/(2.0*epsi)) / (1+sigma(i)*dt/(2.0d0epsi))
    !CHZXL(i) = -(dt/myu) / (1+sigma(i)*dt/(2.0d0*epsi))/dx 
    !CHZY(j)  = (1-sigma(j)*dt/(2.0*epsi)) / (1+sigma(j)*dt/(2.0d0epsi))
    !CHZYL(j) = -(dt/myu) / (1+sigma(j)*dt/(2.0d0*epsi))/dx 

 !    do L=1,6
  !    I0 = LPMLII(L,1)
   !   I1 = LPMLII(L,2)
    !  J0 = LPMLJJ(L,1)
     ! J1 = LPMLJJ(L,2)
 !     K0 = LPMLKK(L,1)
  !    K1 = LPMLKK(L,2)

   ! L1=LPMLST(L)!index要確認
    !do i=I0,I1-1
     !   do j=J0+1,J1-1
      !      do k=K0+1,K1-1
       !         Hxy(L1) = CHXY(i)*Hzx(L1)&
        !                    +CHXYL*(Ez(i,j+1,k)-Ez(i,j,k))
         !       Hxz(L1) = CHXZ(j)*Hzy(L1)&
          !                  +CHXZL*(Ey(i,j,k+1)-Ey(i,j,k)) !要確認
           !     Hx(i,j,k) = Hxy(L1)+Hxz(L1)
            !    L1=L1+1
 !           enddo
  !      enddo
   ! enddo
             
  !  L2=LPMLST(L)
   ! do i=I0+1,I1-1
    !    do j=J0,J1-1
     !       do k=K0+1,K1-1
      !          Hyz(L2) = CHYZ(i)*Hyz(L2)&
       !                     +CHYZL*(Ex(i,j,k+1)-Ex(i,j,k))
        !        Hyx(L2) = CHYX(j)*Hyx(L2)&
         !                   +CHYXL*(Ez(i+1,j,k)-Ez(i,j,k)) !要確認
          !      Hy(i,j,k) = Hyz(L2)+Hyx(L2)
           !     L2=L2+1
  !          enddo
   !     enddo
   ! enddo

   !L3=LPMLST(L)
   ! do i=I0+1,I1-1
    !    do j=J0+1,J1-1
     !       do k=K0+1,K1-1
      !          Hzx(L3) = CHZX(i)*Hzx(L3)&
       !                     +CHZXL*(Ey(i+1,j,k)-Ey(i,j,k))
        !        Hzy(L3) = CHZY(j)*Hzy(L3)&
         !                   +CHZYL*(Ex(i,j+1,k)-Ex(i,j,k)) !要確認
          !      Hz(i,j,k) = Hzx(L3)+Hzy(L3)
           !     L3=L3+1
   !         enddo
    !      enddo
     !  enddo

   ! enddo
    !    endsubroutine H_BoundaryCondition




















































!FFT_fftw3*****************************************************
subroutine dft
    use const_para
    implicit none

    integer :: n,nd,ios
    real(8), allocatable :: acc(:),f(:),p(:)
    complex(kind(0d0)), allocatable :: c(:)
    character(7),parameter :: inp='acc.dat',out1='fas.out',out2='fps.out'
    integer :: plan(8)

    ! FFTW3を呼び出すのに必要なヘッダーファイルを include する
    include 'fftw3.f'

    !ファイル（データ）の長さNDを調べる
    open(51,file=inp,action='read')
    nd=0
    do
        read(51,'(f12.0)',iostat=ios)
        if (ios<0) exit !ファイルの末尾にきたらループを抜ける
         nd=nd+1
    enddo
    close(51)

    !データの長さをnの２乗になるように決めて
    !加速度データ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
    n=2**int(log(dble(nd))/ log(2.0d0)+0.5d0)
    allocate(acc(1:n),c(0:n/2),f(0:n/2),p(0:n/2))

    open(51,file=inp,action='read')
    read(51,'(f12.3)') acc(1:nd)
    close(51)

    !不足分は0パッディング
    acc(nd+1:n)=0.0d0

    !fftwの実行
    call dfftw_plan_dft_r2c_1d(plan,n,acc,c,fftw_estimate)
    call dfftw_execute(plan)
    call dfftw_destroy_plan(plan)

    ! フーリエ振幅スペクトル，フーリエ位相スペクトルを求める
    f= abs(c)
    p= atan2(aimag(c),dble(c))

    ! フーリエ振幅スペクトル，フーリエ位相スペクトルの書き出し
    open(61,file=out1,status='replace',action='write')
    write(61,'(e12.4)') f
    close(61)

    open(61,file=out2,status='replace',action='write')
    write(61,'(f12.5)') p
    close(61)
        endsubroutine dft
























!fftw3_idft***********************************************************************
subroutine idft
    implicit none
    integer :: n,n1,i,ios  !n1いらね iいらね
    real(8), allocatable :: acc(:),f(:),p(:)
    complex(kind(0d0)),allocatable :: c(:),cacc(:)
    character(7),parameter :: inp1='fas.out',inp2='fps.out',out='acc.out'
    integer :: plan(8)

    !fftw3用のヘッダファイルをinclude
    include 'fftw3.f'

    !データの長さNを調べる
    open(51,file=inp1,action='read')
    n=0
    do
     read(51,'(f12.0)',iostat=ios)
     if (ios<0) exit!ファイルの末尾でループをぬける
     n=n+1
    enddo
    close(51)

    !フーリエ逆変換により作成される時刻歴データの長さNを求め、
    !加速度データ、複素フーリエ係数、フーリエ振幅、フーリエ位相の配列を確保
    n=(n-1)*2
    allocate(acc(1:n),cacc(1:n),c(0:n-1),f(0:n/2),p(0:n/2))

    !フーリエ振幅スペクトル、フーリエ位相スペクトルを読み込む
    open(51,file=inp1,action='read')
    read(51,*) f(0:n/2)
    close(51)

    open(51,file=inp2,action='read')
    read(51,*) p(0:n/2)
    close(51)

    !複素フーリエ係数を求める
    c=f*cmplx(cos(p),sin(p))
    c(n/2+1:n-1)=conjg(c(n/2-1:1:-1))

    !fftw3開始
    !planの作成
    call dfftw_plan_dft_1d(plan,n,c,cacc,fftw_backward,fftw_estimate)
    !FFTの実行
    call dfftw_execute(plan)
    !planの破棄
    call dfftw_destroy_plan(plan)

    !フーリエ変換で得られたものを規格化し、実部を取り出す
    acc=dble(cacc/dble(n))

    !加速度時刻歴の書き出し
    open(61,file=out,status='replace',action='write')
    write(61,'(e12.4)') acc
    close(61)
    endsubroutine idft
















!!!FFT 高速フーリエ変換*********************************************
subroutine fft(lx,cx,signi)
    use const_para
    implicit none

    integer, intent(in) :: lx
    integer, intent(in) :: signi
    integer l,m,istep
    real(8) sc
    complex(kind(0d0)),intent(inout) :: cx(lx)
    complex(kind(0d0))  carg,cw,ctemp

   j=1
    sc=sqrt(1./lx)
    do 30 i=1,lx
       if(i>j) go to 10
          ctemp=cx(j)*sc
          cx(j)=cx(i)*sc
          cx(i)=ctemp
10         m=lx/2
20          if(j<=m) go to 30
            j=j-m
          m=m/2
          if(m>=1) go to 20
30             j=j+m
             l=1
40             istep=2*l
             do 50 m=1,l
                carg=(0.,1.)*(3.14159265*signi*(m-1))/l
                cw=cdexp(carg)
                do 50 i=m,lx,istep
                   ctemp=cw*cx(i+l)
                   cx(i+l)=cx(i)-ctemp
50                   cx(i)=cx(i)+ctemp
                   l=istep
                   if(l<lx) goto 40

                     return
                      end subroutine fft




















!!!入力波源の設定  ガウスパルス************************************************
!subroutine gaussianpulse(istep,t,Ie,Mh)
!    use const_para
!    implicit none

!    integer, intent(in) :: istep !タイムステップ
!    real(8), intent(in) :: t     !経過時間
!    real(8), intent(out) :: Ie(nstep) !電流源
!    real(8), intent(out) :: Mh(nstep) !磁流源
!    real(8), parameter :: alpha = (4.0d0/tau0)**2
  
    !0<t<2τ0のときパルス
!    if(t<=2.0d0*tau0) then
!        Ie(istep) = exp(-alpha*(t-tau0)**2)
!        Mh(istep) = exp(-alpha*(t-tau0)**2)

!    else
!        Ie(istep) = 0.0d0
!        Mh(istep) = 0.0d0
!    endif
!            endsubroutine gaussianpulse

