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
!
! 組み込み手続きは必ず"総称名"でつかうこと!
!
!***************************************************

!ficticious wave domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
    use const_para
    implicit none

    integer            :: istep !タイムステップ
    real(8)            :: t !経過時間
    complex(kind(0d0)) :: Je(nstep) !電流源
    complex(kind(0d0)) :: Jh(nstep) !磁流源
!     real(8)            :: sig(nx,ny,nz),myu(nx,ny,nz)
    complex(kind(0d0)) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
!     complex(kind(0d0)) :: Ex(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Ey(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Ez(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hx(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hy(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)) :: Hz(-1:nx+2,-1:ny+2,-1:nz+2)

    open(31,file='hz1000.d')
    open(32,file='hz1010.d')
    open(33,file='hz1020.d')
    open(34,file='hz1030.d')
    open(35,file='hz1040.d')
    open(36,file='hz1050.d')

    open(17,file='je_fic.d')
    open(18,file='jh_fic.d')

    open(20,file='ex1000.d')
    open(21,file='ex1010.d')
    open(22,file='ex1020.d')
    open(23,file='ex1030.d')


open(24,file='hzleft1.d')   !!!反射波の確認
open(25,file='hzright1.d')  !!!反射波の確認
open(26,file='hzleft2.d')
open(27,file='hzright2.d')
open(28,file='hzleft3.d')
open(29,file='hzright3.d')


    !set eh-field to 0
    call set_zero_eh(EX,EY,EZ,HX,HY,HZ)

    write(*,*) '!!!!!!!!!!!!!  start calculation  !!!!!!!!!!!!!!!!'

t=0.0d0!開始時間-----------------------------------

    !モデルの読み込み
    call model
    call init_cpml !cpmlの係数,fdtd部分の係数
    call media_coeff !伝播の係数  　　??cpml以外の媒質のパラメータはどこで設定？

    !cmax,cminの計算 dt,dx,dy,dzの設定
    call confirm_parameter

do istep = 1, nstep !反復計算開始----------------------
!     write(*,*) istep

    !入力波源の設定
    call read_source_3d(istep,t,Hz,Je,Jh)
!     call firstderiv_gauss(istep,t,Je,Jh,sig,myu)
                    !     call read_source_3d(istep,t,sig,myu,Hz,Je,Jh)
!     call read_source_3d(istep,t,sig,myu,EX,Je,Jh)

    !電場計算 E
    call e_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz)
!     call Efield4(istep,t,Je,Ex,Ey,Ez,Hx,Hy,Hz,sig)
!     call Efield4(istep,t,Ex,Ey,Ez,Hx,Hy,Hz,sig)

    !境界条件 CPML_E
    call CPML_E4(ex,ey,ez,hx,hy,hz)!,sig)!
!     call CPML_E2(ex,ey,ez,hx,hy,hz,sig)!
!     call cerjan_e(ex,ey,ez) !!!だめ
!     call cerjan2_e(ex,ey,ez)  !cerjanの吸収境界!!!ok
!     call mur_yz(ex,ey,ez) !Murの吸収境界
!     call mur_zx(ex,ey,ez)!Murの吸収境界
!     call mur_xy(ex,ey,ez)!Murの吸収境界
!     call e_pml(ex,ey,ez,hx,hy,hz) !Normal PML

t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算 H
    call h_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz)
!     call Hfield4(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz,myu)
!     call Hfield4(istep,t,Ex,Ey,Ez,Hx,Hy,Hz,myu)

    !境界条件 CPML_H
    call CPML_H4(ex,ey,ez,hx,hy,hz)!,sig,myu)!,cmax)
!     call CPML_H2(ex,ey,ez,hx,hy,hz,sig,myu)!,cmax)
!     call cerjan_h(hx,hy,hz) !!!だめ
!     call cerjan2_h(hx,hy,hz) !!!ok
!     call h_pml(ex,ey,ez,hx,hy,hz)!Normal PML

t = t + dt*0.5d0 !時間の更新---------------------------

    !アウトプットE-field、H-field
    call output_EH_J(istep,t,Je,Jh,Ex,Ey,Ez,Hx,Hy,Hz)


enddo !*反復計算終了

    !グリーン関数の導出
    !call green()


    close(31)
    close(32)
    close(33)
    close(34)
    close(35)
    close(36)


    close(17)
    close(18)

    close(20)
    close(21)
    close(22)
    close(23)

close(24)
close(25)
close(26)
close(27)
close(28)
close(29)
            end program main
