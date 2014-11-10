fic_dif_mat.f90

module init
  implicit none
  integer,parameter :: s=4096
  integer :: idm
  real(8),parameter ::  pi=3.14159265358979323846d0
  real(8) :: r2,i2,rdm,d
  real(8) :: re2(s),im2(s)
  complex(kind(0d0)):: A(s,s)
  complex(kind(0d0)) :: f_t(s),df(s)
  complex(kind(0d0)),parameter :: complex_j=(0.0d0,1.0d0)
  real(8) :: t
end module


program fic_difmatrix
  use init
  implicit none
  integer :: j,k
  integer :: n,m
  integer :: ts,te,t_rate,t_max,diff
  call system_clock(ts)


  t= 1.14315354400754132E-003

  open(10,file="fic_difmat2.dat")
  do j=1,s
     do k=1,s
        A(j,k) = exp(-(2.0d0*pi*sqrt((j-1)*s*t)*(k-1)/s))*exp(complex_j*(2.0d0*pi*sqrt((j-1)*s*t)*(k-1)/s))
!        A(j,k) =exp(-(2.0d0*pi*sqrt((j)*s*t)*k/s))*exp(complex_j*(2.0d0*pi*sqrt((j)*s*t)*k/s))
        write(10,*) A(j,k)
     end do
     print *, j
  end do
  close(10)


 call system_clock(te,t_rate,t_max)
  if (te<ts) then
     diff=t_max-ts+te
  else
     diff=te-ts
  end if
  print "(A,F10.3)","time it took was:", diff/dble(t_rate)

end program