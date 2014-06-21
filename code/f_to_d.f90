subroutine f_to_d
	use const_para
	implicit none
		do k=0,s-1
			do j=0,s-1


end subroutine f_to_d

!!by kusu
subroutine f_to_d_matrix
	use const_para
	implicit none

		do j=1,s
			do k=1,s
				A(j,k) = exp(-(2.0d0*pi*sqrt((j-1)*s*t)*(k-1/s)))*exp(I_u*(2.0d0*pi*sqrt((j-1)*s*t)*(k-1)/s))
			enddo
		enddo
end subroutine f_to_d_matrix

!by imamu
subroutine laplace_fft
	use const_para
	implicit none
		integer :: n
		real(8) :: om
		om =2.0d0*pi/it/dt
		t0=pi/fmax_w
		beta=pi*fmax**2.0d0


	do n=1,it

		do k=1,it
			EX_w(n) = EX_w(n) &
					 + EX_f(k)*dt *exp(-sqrt(omega0*om*n)*k*dt) *cexp(I_u*sqrt(omega0*om*n)*k*dt)

			JX_w(n) = JX_w(n) &
					+ csqrt(-2.0d0*omega0/I_u/om/n)*JX_f(k)*dt *exp(-sqrt(omega0*om*n)*k*dt) *cexp(I_u*sqrt(omega0*om*n)*k*dt)
		enddo
			JX_w(0) = 2.0d0 * omega0  !!!要確認

			GX_w(n) = EX_w(n) / JX_w(n)

	enddo
end subroutine laplace_fft



subroutine convolution_GJ_to_E
	use const_para
	implicit none
		complex(kind(0d0)) :: in_G
		complex(kind(0d0)) :: in_J
		complex(kind(0d0)) :: in_EF
		complex(kind(0d0)) :: out_G
		complex(kind(0d0)) :: out_J
		complex(kind(0d0)) :: out_ET

end subroutine convolution_GJ_to_E














