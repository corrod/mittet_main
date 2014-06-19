
program sin

  implicit none
  integer :: i,n
  real(8) :: func_sin
  real(8) :: dx,dy,dz,dt,t,t0,freq,la,mu0,eps0,pi,c,tp,t1
  c=3.0d8
  dx=1.0d-2
  dy=dx
  dz=dx
  freq=2.45d9
!  t0=0.646d0/freq
!  la=(1.0d0/0.29d0/t0)**2
  n=1000
  pi=3.141592
  mu0=4.0d-7*pi
  eps0=8.854d-12
  tp=1.0d0/freq
  t1=3.0d0*tp

  dt=1.0d0/freq/1.0d2
  do i=1,n
     t=0.0d0+dble(i-1)*dt
 ! if(t.lt.tp) then
     !     : exp(-(x-m)^2/(2 s^2))

     if(t.lt.t1) then
        func_sin=0.5d0*(1.0d0-dcos(pi*t/t1))*dsin(2.0d0*pi*freq*t)
     else

        func_sin=dsin(2.0d0*pi*freq*t)
  end if
  open(30,file='sin.dat')
  write(30,*) i*dt,  func_sin
  end do
  close(30)
end program
