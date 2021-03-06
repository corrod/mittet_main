!複素２次元離散フーリエ変換

!実数データをFFT
program 2dfft_r2c
    use const_para
    implicit none

    integer, parameter :: M=32, N=32
    real(8)  :: in(M,N)
    complex(kind(0d0)) :: out(M/2+1,N)
    integer(8) :: plan

    include 'fftw3.f'


    call dfftw_plan_dft_r2c_2d(plan, M, N, in, out, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, in, out)


    open(11,file='out.dat',form='formatted')
        do i=1,N
            write(11,*) i, real(out(i)), aimag(out(i))
        enddo
    close(11)

    call dfftw_destroy_plan(plan)
end program 2dfft_r2c


!次元
program fft_2d_fw
    use const_para
    implicit none

    integer, parameter :: M=32, N=32
    complex(kind(0d0)) :: in(M,N)
    complex(kind(0d0)) :: out(M,N)
    integer(8) :: plan

    include 'fftw3.f'

    call dfftw_plan_dft_2d(plan, M, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, in, out)

    open(11,file='out.dat',form='formatted')
        do i=1,N
            write(11,*) i, real(out(i)), aimag(out(i))
        enddo
    close(11)

    call dfftw_destroy_plan(plan)
end program fft_2d_fw

!２次元
program fft_2d_bw
    use const_para
    implicit none

    integer, parameter :: M=32, N=32
    complex(kind(0d0)) :: in(M,N)
    complex(kind(0d0)) :: out(M,N)
    integer(8) :: plan

    include 'fftw3.f'

    call dfftw_plan_dft_2d(plan, M, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, in, out)

    open(11,file='out.dat',form='formatted')
        do i=1,N
            write(11,*) i, real(out(i)), aimag(out(i))
        enddo
    close(11)

    call dfftw_destroy_plan(plan)
end program fft_2d_fw
