
program rironkai

  implicit none
  integer :: i,n
  real(8) :: theory
  real(8) :: l,dt,t,t0,freq,la,mu0,eps0,pi,c,tp,t1,r,j
  j=1.0d-4
  c=3.0d8
  l=1.0d-2
  freq=2.45d9
  n=1000
  pi=3.141592
  mu0=4.0d-7*pi
  eps0=8.854d-12
  tp=1.0d0/freq
  t1=3.0d0*tp
  r=4.8d-1

  dt=1.0d0/freq/1.0d2
  do i=1,n
     t=0.0d0+dble(i-1)*dt-r/c

     if(t.lt.0) then
     theory=0.0d0
     else if(t.lt.t1) then
        theory=0.5d0*(1.0d0-dcos(pi*t/t1))*(6.0d1*pi*l*j)/((c/freq)*r)*dabs(dsin(2.0d0*pi*freq*t-pi/2.0d0))
     else

        theory=(6.0d1*pi*l*j)/((c/freq)*r)*dabs(dsin(2.0d0*pi*freq*t-pi/2.0d0))
  end if
  open(30,file='theory.dat')
  write(30,*) i*dt,  theory
  end do
  close(30)
end program
