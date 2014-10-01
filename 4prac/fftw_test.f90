program main
    implicit none
	integer :: i
	integer(8) :: x(64)
	complex(kind(0d0)) :: in1(64)
!real(8) :: in1(64)
	complex(kind(0d0)) :: out1(64)
	integer :: plan(8)

include 'fftw3.f'


open(44,file='ffttestinput.d')
do i=1,64
read(44,*) x(i), in1(i)
enddo
close(44)

do i=1,64
write(10,*) i-1,real(in1(i)),aimag(in1(i))
!write(10,*) i-1, in1(i)
enddo

!////////////////////////////////////////////////////////////
! make plans
!      FFTW_FORWARD (-1) or FFTW_BACKWARD (+1)
!////////////////////////////////////////////////////////////
	call dfftw_plan_dft_1d(plan,64,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE) !complex array入力
!	call dfftw_plan_dft_1d(plan,64,in1,out1,FFTW_BACKWARD,FFTW_ESTIMATE) !complex array入力



!///////////////////////////////////////////////////////////
! carry out fourier transformation
!///////////////////////////////////////////////////////////
	call dfftw_execute(plan,in1,out1)
!	call dfftw_execute(plan)


!///////////////////////////////////////////////////////////
! destroy plan
!///////////////////////////////////////////////////////////
	call dfftw_destroy_plan(plan)
do i=1,64
out1(i)=out1(i) / 64.0d0
write(11,*) i,real(out1(i)),aimag(out1(i))
enddo

do i=1,64
write(12,*) i,abs(out1(i))
enddo

end program main

