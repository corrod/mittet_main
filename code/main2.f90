!***************************************************
!  mittet 仮想領域法 2014.06.01
!
!
!   omega0=2πf0
!   f0=1,0Hz
!   sigwa=3.2S/m
!
!   x=(i-1)*dx
!   代入必要?↓↓↓
!   E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!***************************************************


!ficticious wave domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
    use const_para
    implicit none

    integer            :: istep !タイムステップ
    real(8)            :: t !経過時間
    real(8)            :: Je(nstep) !電流源
    real(8)            :: Jh(nstep) !磁流源
    real(8)            :: sig(nx,ny,nz),myu(nx,ny,nz)
    complex(kind(0d0)) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
!     complex(kind(0d0)) :: Ex(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Ey(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Ez(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hx(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hy(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hz(-1:nx+2,-1:ny+2,-1:nz+2)

    open(13,file='hz1100.d')
    open(14,file='hz1010.d')
    open(15,file='hz1020.d')
    open(16,file='hz1030.d')
    open(17,file='jh_fic.d')

write(*,*) '!!!!!!!!!!!!!  start calculation  !!!!!!!!!!!!!!!!'

    t=0.0d0!開始時間---------------------------------------

!   Ex,Ey,Ez,Hx,Hy,Hzの初期化
!     Ex(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
!     Ey(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
!     Ez(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
!     Hx(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
!     Hy(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
!     Hz(-1:nx+2,-1:ny+2,-1:nz+2)=0.0d0
    Ex(1:nx,1:ny,1:nz)=0.0d0
    Ey(1:nx,1:ny,1:nz)=0.0d0
    Ez(1:nx,1:ny,1:nz)=0.0d0
    Hx(1:nx,1:ny,1:nz)=0.0d0
    Hy(1:nx,1:ny,1:nz)=0.0d0
    Hz(1:nx,1:ny,1:nz)=0.0d0


    !モデルの読み込み
    call model(sig,myu)

    !cmax,cminの計算 dt,dx,dy,dzの設定
    call set_d_txyz !(cmax)

do istep = 1, nstep !反復計算開始----------------------
write(*,*) istep

    !入力波源の設定
    call firstderiv_gauss(istep,t,Je,Jh,sig,myu)

    !電場計算 E
    call Efield(istep,t,Je,Ex,Ey,EZ,Hx,Hy,Hz,sig)

    !境界条件 CPML_E
    call CPML_E(ex,ey,ez,hx,hy,hz,sig)!,cmax)  !
!     call cerjan_e(ex,ey,ez) !!!だめ
!     call cerjan2_e(ex,ey,ez)  !cerjanの吸収境界!!!ok
!     call mur_yz(ex,ey,ez) !Murの吸収境界
!     call mur_zx(ex,ey,ez)!Murの吸収境界
!     call mur_xy(ex,ey,ez)!Murの吸収境界
!     call e_pml(ex,ey,ez,hx,hy,hz) !Normal PML

t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算 H
    call Hfield(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz,myu)

    !境界条件 CPML_H
    call CPML_H(ex,ey,ez,hx,hy,hz,sig,myu)!,cmax)
!     call cerjan_h(hx,hy,hz) !!!だめ
!     call cerjan2_h(hx,hy,hz) !!!ok
!     call h_pml(ex,ey,ez,hx,hy,hz)!Normal PML

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
