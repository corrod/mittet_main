
program hagen

  implicit none
  integer :: i,n
  real(8) :: func_gaussian1,func_gaussian
  real(8) :: dx,dy,dz,dt,t,t0,freq,la,mu0,eps0,pi,c,tp,t1,beta,omg
  real(8) :: eps_p(3),c_p(3),sig(3)
  c=3.0d8
  dx=1.0d2
  dy=dx
  dz=dx
  freq=1.0d0

  n=200
  pi=3.141592
  mu0=4.0d-7*pi
  eps0=8.854d-12
  sig(3)=3.2d0


  omg=2.0d0*pi*freq
  eps_p(3)=sig(3)/(2.0d0*omg)
  c_p(3)=sqrt((2.0d0*omg)/(mu0*sig(3)))

  tp=pi/freq
  beta=pi*freq**2
  dt=0.99*(1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(3))
 
  do i=1,n
     t=0.0d0+dble(i-1)*dt

!  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)
  func_gaussian=dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)

  open(30,file='hagen.dat')
  write(30,*) i,  func_gaussian
  end do
  close(30)
end program
