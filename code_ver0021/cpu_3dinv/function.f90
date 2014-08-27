! 各種関数 function

!残差計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine residualE()
    implicit none

    character :: res1(256), res2(256), res3(256)
    complex(kind(0d0)) :: delEx_f(), delEy_f(), delEz_f()
    complex(kind(0d0)) :: Exobs(), Eyobs(), Ezobs()
    complex(kind(0d0)) :: EcalX_ff(), EcalY_ff(), EcalZ_ff()
    real(8) :: err_sum(MAXITER)
    !出力ファイルの作成
    open(1,file='residualEX.dat')
    open(2,file='residualEY.dat')
    open(3,file='residualEX.dat')

    !残差計算
    delEx_f() = Exobs() - EcalX_ff()
    delEy_f() = Eyobs() - EcalY_ff()
    delEz_f() = Ezobs() - EcalZ_ff()
    !出力
    write(1,*)
    write(2,*)
    write(3,*)

    !error 合計計算
    err_sum(iter) = err_sum(iter) + abs(delEx_f(isr)) !itr = itr回数
    err_sum(iter) = err_sum(iter) + abs(delEy_f(isr)) !  isr  = irec*upperF + w
    err_sum(iter) = err_sum(iter) + abs(delEz_f(isr))

    close(1)
    close(2)
    close(3)

end subroutine residualE





!sensitivityの計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine sensitivity()
    implicit none

    sensitiv = sensitiv + abs(EcalX_ff * EcalX_bb)
    sensitiv = sensitiv + abs(EcalY_ff * EcalY_bb)
    sensitiv = sensitiv + abs(EcalZ_ff * EcalZ_bb)
end subroutine sensitiviry





!  2.6 delE を green fucntion と convolutoin する
!     EcalX_back = EcalX_back + EcalX_bb * conjg(delEx_f)

!delE を green fucntion と convolutoin する!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine conv_GdelE()
    implicit none

    EcalX_back() = 0.0d0
    EcalY_back() = 0.0d0
    EcalZ_back() = 0.0d0




    EcalX_back() = EcalX_back() + EcalX_bb() * conjg(delEx_f())
    EcalY_back() = EcalY_back() + EcalY_bb() * conjg(delEx_f())
    EcalZ_back() = EcalZ_back() + EcalZ_bb() * conjg(delEx_f())

    EcalX_back() = EcalX_back() + EcalX_bb() * conjg(delEy_f())
    EcalY_back() = EcalY_back() + EcalY_bb() * conjg(delEy_f())
    EcalZ_back() = EcalZ_back() + EcalZ_bb() * conjg(delEy_f())

    EcalX_back() = EcalX_back() + EcalX_bb() * conjg(delEz_f())
    EcalY_back() = EcalY_back() + EcalY_bb() * conjg(delEz_f())
    EcalZ_back() = EcalZ_back() + EcalZ_bb() * conjg(delEz_f())

end subroutine conv_GdelE





!gammaの計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine calc_gamma()
    implicit none

    open(3,file='gradient.dat')

!do loop
!gradient の計算
    grad() = grad() + real(EcalX_ff() * EcalX_back() &
                     + EcalY_ff() * EcalY_back() &
                     + RcalZ_ff() * EcalZ_back()) / sensitiv()
    write(3,*) i*dx, j*dy, k*dz, grad()

    close(3)

end subroutine calc_gamma




!calc phi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine calc_phi()
    implicit none

    real(8) :: sigmax2
    real(8) :: gradmax
    real(8) :: sig2

    do k=1,nz
        do j=1,ny
            do i=1,nx
                if(k<=isebed) then !　
                    if(sig() > sigmax2) sigmax2 = sig() !　
                    if(grad() > gradmax) gradmax = grad() !　
                endif
            endif
        endif
            enddo
        enddo
    enddo

    scaler = 0.01 * sigmax2 / gradmax

    open(1,file='sig_tmp.dat')

    do k=1,nz
        do j=1,ny
            do i=1,nx
                sig2() = sig() - scaler * grad()
                write(1,*) i*dx, j*dy, k*dz, sig2()
            enddo
        enddo
    enddo


    close(1)

end subroutine calc_phi



!calc alpha!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine calc_alpha()
    implicit none

    complex(kind(0d0)) :: ctmp_x, ctmp_y, ctmp_z
    complex(kind(0d0)) :: EcalX_ff2(), EcalY_ff2(), EcalZ_ff2()
    complex(kind(0d0)) :: EcalX_ffL(), EcalY_ffL(), EcalZ_ffL()
    complex(kind(0d0)) :: k1, k2
!doloop
    ctmp_x = EcalX_ff2() - EcalX_ffL()
    ctmp_y = EcalY_ff2() - EcalY_ffL()
    ctmp_z = EcalZ_ff2() - EcalZ_ffL()

    k1 = k1 + delEx_f() * ctmp_x
    k1 = k1 + delEy_f() * ctmp_y
    k1 = k1 + delEz_f() * ctmp_z

    k2 = k2 + ctmp_x * ctmp_x
    k2 = k2 + ctmp_y * ctmp_y
    k2 = k2 + ctmp_z * ctmp_z

end subroutine calc_alpha





!update parameter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine update_para()
    implicit none

    open(1, file = "gradient_03d.dat")
    open(2, file = "sig_03d.dat")

    alpha = 3.0d0 * abs(k1) / abs(k2)  !why 3 ?

!doloop
    if(k<=iseabed) then
    sig() = sig() + alpha * scaler * grad()
    elseif(sig()<minval) sig() = minval
    elseif(sig()>maxval) sig() = maxval

!doloop
    do k=1,nz
        do j=1,ny
            do i=1,nx
                write(1,*) i*dx, j*dy, k*dz, alpha*scaler*grad()
                write(2,*) i*dx, j*dy, k*dz, sig()
                grad_p() = grad()
                grad() = 0.0d0
            enddo
        enddo
    enddo

    close(1)
    close(2)

end subroutine update_para


!media coeff sig tmp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine media_coeff_sig_tmp()
    implicit none

    real(8) :: eps2

    open(1,file='coef1.dat')
    open(2,file='coef2.dat')
    open(3,file='coef3.dat')

    do k=1,nz
        do j=1,ny
            do i=1,nx

                eps2 = sig2() / 2.0d0 / omega0

                !CPML coeff1
                ca_x(i,j,k) = (1.0d0 - ((esigx(i)*dt)/(2.0d0*eps2))) &
                            / (1.0d0 + ((esigx(i)*dt)/(2.0d0*eps2)))
                ca_y(i,j,k) = (1.0d0 - ((esigy(j)*dt)/(2.0d0*eps2))) &
                            / (1.0d0 + ((esigy(j)*dt)/(2.0d0*eps2)))
                ca_z(i,j,k) = (1.0d0 - ((esigz(k)*dt)/(2.0d0*eps2))) &
                            / (1.0d0 + ((esigz(k)*dt)/(2.0d0*eps2)))
                !CPML coeff2
                da_x(i,j,k) = (1.0d0 - ((msigx(i)*dt)/(2.0d0*eps2))) $
                            / (1.0d0 + ((msigx(i)*dt)/(2.0d0*eps2)))
                da_y(i,j,k) = (1.0d0 - ((msigy(j)*dt)/(2.0d0*eps2))) $
                            / (1.0d0 + ((msigy(j)*dt)/(2.0d0*eps2)))
                da_z(i,j,k) = (1.0d0 - ((msigz(k)*dt)/(2.0d0*eps2))) $
                            / (1.0d0 + ((msigz(k)*dt)/(2.0d0*eps2)))
                !CPML coeff3
                cb_x(i,j,k) = dt / eps2 / (1.0d0+(esigx(i)*dt)/(2.0d0*eps2))
                cb_y(i,j,k) = dt / eps2 / (1.0d0+(esigy(j)*dt)/(2.0d0*eps2))
                cb_z(i,j,k) = dt / eps2 / (1.0d0+(esigz(k)*dt)/(2.0d0*eps2))
                !CPML coeff4
                db_x(i,j,k) = dt / MU0 / (1.0d0 + (msigx(i)*dt)/(2.0d0*eps2))
                db_y(i,j,k) = dt / MU0 / (1.0d0 + (msigy(j)*dt)/(2.0d0*eps2))
                db_z(i,j,k) = dt / MU0 / (1.0d0 + (msigz(k)*dt)/(2.0d0*eps2))

                !　　　
                write(1,*) ca_x(i,j,k), ca_y(i,j,k), ca_z(i,j,k),cb_x(i,j,k),cb_y(i,j,k),cb_z(i,j,k)
                write(2,*) da_x(i,j,k), da_y(i,j,k), da_z(i,j,k),db_x(i,j,k),db_y(i,j,k),db_z(i,j,k)
                write(3,*) esigx(i),esigy(j),esigz(k),msigx(i),msigy(j),msigz(k)
            enddo
        enddo
    enddo

    close(1)
    close(2)
    close(3)

end subroutine media_coeff_sig_tmp




!laplace to freq!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!　
subroutine laplaceToFreq()
    implicit none

    real(8) :: om, t0, beta
    complex(kind(0d0)) :: ctmp_f1(), ctmp_f2(), ctmp_f3()
    complex(kind(0d0)) :: ctmp_ex(), ctmp_ey, ctmp_ez
    complex(kind(0d0)) :: ctmp_j
    complex(kind(0d0)) :: ctmp_gx, ctmp_gy, ctmp_gz
    integer(8) :: plan1,

include fftw3 !　　　　　　　　　　　　　　　fftw3をついが

    om = 2.0d0 * M_PI / it / dt
    t0 = M_PI / fmax_w
    beta = M_PI * (fmax_w**2.0d0)


!　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
call dfftw_plan_dft_1d(plan1,nd,in1,out1,FFTW_BACKWARD,fftw_estimate)
call dfftw_execute(plan1,in1,out1)
call dfftw_destroy_plan(plan1)


!doloop
    ctmp_f1() = EcalX_b()
    ctmp_f2() = EcalY_b()
    ctmp_f3() = EcalZ_b()


    !laplace transform
    do w =1, !　
        ctmp_ex(w) = 0d0
        ctmp_ey(w) = 0d0
        ctmp_ez(w) = 0d0

        ctmmp_j(w) = 0d0
    do k=1,it
        ctmp_ex(w) = ctmp_ex(w) + ctmp_f1(k) * dt &
                     * exp(-sqrt(omega0*om*w)*k*dt) * exp(I_u*sqrt(omega0*om*w)*k*dt)
        ctmp_ey(w) = ctmp_ey(w) + ctmp_f2(k) * dt &
                     * exp(-sqrt(omega0*om*w)*k*dt) * exp(I_u*sqrt(omega0*om*w)*k*dt)
        ctmp_ez(w) = ctmp_ez(w) + ctmp_f3(k) * dt &
                     * exp(-sqrt(omega0*om*w)*k*dt) * exp(I_u*sqrt(omega0*om*w)*k*dt)
        ctmp_j(w) = ctmp_j(w) + sqrt(-2.0d0*omega0/I_u/om/w) * JX_f(k) * dt &
                     * exp(-sqrt(omega0*om*w)*k*dt) * exp(I_u*sqrt(omega0*om*w)*k*dt)

        ctmp_j(0) = 2.0d0 * omega0

        ctmp_gx(w) = ctmp_ex(w) / ctmp_j(w)
        ctmp_gy(w) = ctmp_ey(w) / ctmp_j(w)
        ctmp_gz(w) = ctmp_ez(w) / ctmp_j(w)
        !convolution G * J
        EcalX_ff2() = ctmp_gx(w) * out_True(w)
        Ecaly_ff2() = ctmp_gy(w) * out_True(w)
        EcalZ_ff2() = ctmp_gz(w) * out_True(w)

end subroutine laplaceToFreq
