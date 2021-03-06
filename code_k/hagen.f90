
program hagen

  implicit none
  integer :: i,n
  real(8) :: func_gaussian1
  real(8) :: dx,dy,dz,dt,t,t0,freq,la,mu0,eps0,pi,c,tp,t1,beta
  c=3.0d8
  dx=1.0d-2
  dy=dx
  dz=dx
  freq=2.45d9
!  t0=0.646d0/freq
!  la=(1.0d0/0.29d0/t0)**2
  n=100
  pi=3.141592
  mu0=4.0d-7*pi
  eps0=8.854d-12
  t1=3.0d0*tp

  tp=pi/freq
  beta=pi*freq**2

  dt=1.0d0/freq/20.0d0
  do i=1,n
     t=0.0d0+dble(i-1)*dt

  func_gaussian1=-2.0d0*beta*(t-tp)*dsqrt(beta/pi)*dexp(-beta*(t-tp)**2)

  open(30,file='hagen.dat')
  write(30,*) i*dt,  func_gaussian1
  end do
  close(30)
end program
