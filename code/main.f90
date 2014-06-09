!***************************************************
!  mittet 仮想領域法 2014.06.01
!
!
!
!   omega0=2πf0
!   f0=1,0Hz
!   sigmawa=3.2S/m
!   x=(i-1)*dx
!   operator half length ∂n=Σαの書き方
!   代入必要?↓↓↓
!   E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!***************************************************

!!!メインプログラム************************************
!****************************************************
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
    open(16,file='hz0.d') !iran
    open(13,file='hz10.d') !iran
    open(14,file='hz20.d') !iran
    open(17,file='jh.d')!iran

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

    do istep = 1, nstep !*反復計算開始----------------------

    !入力波源の設定
!   call gaussianpulse(istep,t,Ie,Mh)
    call gaussian(istep,t,Je,Jh,sigma,myu)
    !電場計算
    call EXFIELD(istep,t,Je,Ex,Hy,Hz,sigma)
    call EYFIELD(istep,t,Je,Ey,Hz,Hx,sigma)
    call EZFIELD(istep,t,Je,Ez,Hx,Hy,sigma)

     !境界条件 E
!   call E_PML(Ex,Ey,Ez,HX,Hy,Hz,sigma)
!    call CPML_E(ex,ey,ez,hx,hy,hz)

    t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算
    call HXFIELD(istep,t,Jh,Hx,Ey,Ez,myu)
    call HYFIELD(istep,t,Jh,Hy,Ex,Ez,myu)
    call HZFIELD(istep,t,Jh,Hz,Ex,Ey,myu)

    !境界条件 H
!   call H_BoundaryCondition(Ex,Ey,Ez,Hx,Hy,Hz,sigma,myu)
 !   call CPML_H(ex,ey,ez,hx,hy,hz)

    t = t + dt*0.5d0 !時間の更新---------------------------

    !アウトプットE-field、H-field
    call output_EH(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz)

    enddo !*反復計算終了


    !グリーン関数の導出
    !call green()

    !高速フーリエ変換
   ! call ficticiou_to_diffusive_f(t)
    close(16)
    close(13)
    close(14)
    close(17)
            end program main
