program makewave
	use const_para
		implicit none
! 		real(8) :: fmax_w, vc
! 		complex(kind(0d0)) :: JX_t,JX_f,JX_w
		real(8) :: dtmp
		character(8) :: invJT,invJF,invJW
! 		invJT=
! 		invJF=
! 		invJW=
		call init_var()
		call init_emfield()

		call fourier_trans(invJT,invJF,invJW)
		call inv_fourier(invJT,invJF,invJW)
end program makewave


subroutine init_var()

		implicit none
		integer :: tmp
		real(8) :: dtmp
		real(8) :: minv,maxv
		real(8) :: cmax,cmin
		real(8) :: courant
		character(8) :: ifm, ifm2

		ifm =
		ifm2=

		cmin = sqrt(2.0d0*omega0/MU0/maxv)
		cmax = sqrt(2.0d0*omega0/MU0/minv)
		courant = 1.0d0/cmax/sqrt(1.0d0/dx**2.0d0 + 1.0d0/dy**2.0d0 + 1.0d0/dz**2.0d0)
		dt = courant*6.0d0/7.0d0*0.999d0
	end subroutine init_var


subroutine init_emfield()
	implicit none
		integer :: i
		real(8) :: dtmp,ReF,ImF
		character(8) :: ifw1,ifw2,ifw3
 		complex(kind(0d0)) :: JX_t,JX_f,JX_w

 		ReF=0.0d0
 		ImF=0.0d0

 		do i=1,N
 			open()
 			read() dtmp, ReF
 			JX_t(i) = ReF
 			close()
 		end subroutine init_emfield



subroutine fourier_trans()
	implicit none
	integer :: n,nstep
	real(8) :: om
	om = 2.0d0*pi/it/dt   !!! it?
	do n=1,N
		om = 2.0d0*pi/N/dt
		Jx_f(n) = 0.0d0
		do tt=1,N
			JX_f(n) = JX_f(n) + JX_f(n) *exp(I_u*2.0d0*pi*tt*n/dble(N))
		enddo
			JX_w(n) = sqrt(-I_u*om*dble(n)/2.0d0/omega0) * JX_f(n)
	enddo

	open(1,file='inv_JF')
		write(1) om*n, real(JX_f(n)),real(I_u*JX_f(n)),real(JX_f(n))*real(JX_f(n))+real(I_u*JX_f(n))*real(I_u*JX_f(n)),atan(real(JX_f(n)*I_u)/real(JX_f(n)))
		close(1)

	open(2)
		write(2,file='inv_JW') om*n, real(JX_w[n]),real(I_u*JX_w[n]), real(JX_w[n])*real(JX_w[n])+real(I_u*JX_w[n])*real(I_u*JX_w[n]),atan(real(JX_w[n]*I_u)/real(JX_w[n]))
		close(2)
	end subroutine fourier_trans


subroutine inv_fourier
	implicit none
	integer :: n,k
	real(8) :: om

	do k=1,N
		om = 2.0d0*pi/N/dt
		do n=1,N
			JX_t(k) =JX_t(k) + (I_u+1) * sqrt(-I_u/2.0d0)/4.0d0/pi+JX_f(n)*exp(sqrt(om*n*omega0)*k*dt)*exp(-I_u*sqrt(om*dble(n)*omega0)*dble(k)*t)

	open()
	write() dt*k,real(JX_t(k)),real(I*JX_t(k))
	close()














