!***************************************************
!  mittet 仮想領域法 2014.06.01
!
!
!   omega0=2πf0
!   f0=1,0Hz
!   sigmawa=3.2S/m
!
!   x=(i-1)*dx
!   代入必要?↓↓↓
!   E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!***************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!      メインプログラム     !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
    use const_para
    implicit none

    integer            :: istep !タイムステップ
    real(8)            :: t !経過時間
    real(8)            :: Jn(nstep) !gaussian
    real(8)            :: Je(nstep) !電流源
    real(8)            :: Jh(nstep) !磁流源
    real(8)            :: sigma(nx,ny,nz),myu(nx,ny,nz)
    real(8)            :: cmax
    complex(kind(0d0)) :: Ex(nx,ny,nz)
    complex(kind(0d0)) :: Ey(nx,ny,nz)
    complex(kind(0d0)) :: Ez(nx,ny,nz)
    complex(kind(0d0)) :: Hx(nx,ny,nz)
    complex(kind(0d0)) :: Hy(nx,ny,nz)
    complex(kind(0d0)) :: Hz(nx,ny,nz)

    open(13,file='hz1000.d') 
    open(14,file='hz1010.d') 
    open(15,file='hz1020.d') 
    open(16,file='hz1030.d') 
    open(17,file='jh_fic.d')
    
write(*,*) '!!!!!!!!!!!!!  start calculation  !!!!!!!!!!!!!!!!'

    t=0.0d0!開始時間---------------------------------------

!   Ex,Ey,Ez,Hx,Hy,Hzの初期化
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                Ex(i,j,k)=0.0d0
                Ey(i,j,k)=0.0d0
                Ez(i,j,k)=0.0d0
                Hx(i,j,k)=0.0d0
                Hy(i,j,k)=0.0d0
                Hz(i,j,k)=0.0d0
            enddo
        enddo
    enddo

    !モデルの読み込み
    call model(sigma,myu)

    !cmax,cminの計算 dt,dx,dy,dzの設定
    call set_d_txyz(cmax)


    do istep = 1, nstep !反復計算開始----------------------

    !入力波源の設定
    call gaussian(istep,t,Je,Jh,sigma,myu)

    !電場計算 E
    call Efield(istep,t,Je,Ex,Ey,EZ,Hx,Hy,Hz,sigma)
    
    !境界条件 CPML_E
    call CPML_E(ex,ey,ez,hx,hy,hz,sigma,cmax)


    t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算 H
    call Hfield(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz,myu)

    !境界条件 CPML_H
    call CPML_H(ex,ey,ez,hx,hy,hz,sigma,myu,cmax)


    t = t + dt*0.5d0 !時間の更新---------------------------

    !アウトプットE-field、H-field
    call output_EH(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz)

    enddo !*反復計算終了

    !グリーン関数の導出
    !call green()


    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
            end program main
