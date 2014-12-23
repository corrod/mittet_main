subroutine residual_bp

end subroutine residual_bp


program main

    use const_para
    implicit none

    integer            :: istep !タイムステップ
    real(8)            :: t !経過時間
    complex(kind(0d0)) :: Je(nstep) !電流源
    complex(kind(0d0)) :: Jh(nstep) !磁流源
    complex(kind(0d0)) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)


do sp = 1+ncpml+L_hori,nx-ncpml-L_hori

    write(file_sp,*) sp !

    open(15,file="signal_sin.d")
    open(16,file='signal.d')
    open(17,file='je_fic.d')
    open(18,file='jh_fic.d')


    open(55,file = trim(adjustl(file_sp))//'pattern2_1.d')   !パターン２
    open(56,file = trim(adjustl(file_sp))//'pattern2_2.d')   !①
    open(57,file = trim(adjustl(file_sp))//'pattern2_3.d')   !② ③


    write(*,*) 'x_source position :', sp

    !set EM-field to 0
    call set_zero_eh(EX,EY,EZ,HX,HY,HZ)

    write(*,*) '!!!!!!!!!!!!!  start calculation  !!!!!!!!!!!!!!!!'

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
t=0.0d0!開始時間-----------------------------------

    !送受信位置設定
    call source_posi

    !モデルの読み込み
!     call model_nocrack !欠陥なし
    call model_simple !1つのピンホール
!     call model_stair !階段状モデル1
!     call model_pinhole8 !深さ8mmの場合のピンホールモデル
!     call model_pinholeall !全部のピンホールモデル

    !cpmlの係数,fdtd部分の係数
    call init_cpml

    !伝播fdtdの係数
    call media_coeff !　??cpml以外の媒質のパラメータはどこで設定？

    !cmax,cminの計算 dt,dx,dy,dzの設定
    call confirm_parameter

do istep = 1, nstep !反復計算開始----------------------

    !入力波源の設定
    call read_source_3d(istep,t,Hz,Je,Jh)

    !電場計算 E
    call e_field_cpml4bp(Ex,Ey,EZ,Hx,Hy,Hz)

    !境界条件 CPML_E
    call CPML_E4(ex,ey,ez,hx,hy,hz)

t = t + dt*0.5d0  !時間の更新--------------------------

    !磁場計算 H
    call h_field_cpml4bp(Ex,Ey,EZ,Hx,Hy,Hz)

    !境界条件 CPML_H
    call CPML_H4(ex,ey,ez,hx,hy,hz)

t = t + dt*0.5d0 !時間の更新---------------------------

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


    !アウトプットE-field、H-field
    call output_EH_J(istep,t,Je,Jh,Ex,Ey,Ez,Hx,Hy,Hz)

enddo !*反復計算終了

    close(15)
    close(16)
    close(17)
    close(18)


    close(55) !diff_probe用
    close(56)
    close(57)

    !差分プローブ 自己比較、相互比較
    call diff_probe

    !ficticious から diffusiveへ E
    call f_to_d_e

    !ficticious から diffusiveへ H
    call f_to_d_h

enddo !ソース位置 loop

        end program main