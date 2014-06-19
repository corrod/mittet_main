
program hagen

  implicit none
  integer :: i,n
  real(8) :: func_gaussian1,func_gaussian,func_sin
  real(8) :: dx,dy,dz,dt,t,t0,freq,la,mu0,eps0,pi,c,tp,t1,beta,omg
  real(8) :: eps_p(7),c_p(7),sig(7)
  c=3.0d8
  dx=2.0d1
  dy=dx
  dz=dx
  freq=1.0d0

  n=16384
  pi=3.1415926539
  mu0=4.0d-7*pi
  eps0=8.854d-12
 ! sig(3)=3.2d0
  sig(7)=1.0d-1

  omg=2.0d0*pi*freq
  eps_p(7)=sig(7)/(2.0d0*omg)
  c_p(7)=sqrt((2.0d0*omg)/(mu0*sig(7)))

  tp=pi/(freq*12.5d0)
  beta=pi*(12.5d0*freq)**2
  dt=0.99*(1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(7))
 
  do i=0,n-1
     t=0.0d0+dble(i)*dt

!  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)/237.26829917425468
  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)


!  func_gaussian=dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)
!     func_sin=dsin(2.0d0*pi*(freq*1.0d1)*t)

  open(30,file='func_gauss1f.dat')
  write(30,*) i,i*dt,  func_gaussian1
  end do
  close(30)
end program
