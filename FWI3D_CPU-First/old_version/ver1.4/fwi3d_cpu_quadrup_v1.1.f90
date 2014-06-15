program quad
  implicit none
  real(16),parameter :: PI =4.d0*atan(1.d0)
  real(16),parameter :: MU0=4.d0*PI*1.e-7
  real(16),parameter :: sx=0.d0, sy=0.d0, sz=0.d0
  real(16),parameter :: rx=5000.d0, ry=0.d0, rz=0.d0
  real(16),parameter :: moment =1.d0
  real(16),parameter :: sig =1.d0

  integer NN,t
  real(16) :: dx,dy,dz, dt, omega_0,fmax_w,vc
  real(16) :: coef1,coef2,time,r,t1,t0,beta
  real(16) :: dtmp, courant,om,ReF
  real(16) :: EX,EX_p,d_EX,gamma1,d_gamma,d2_gamma
  real(16) :: theta, theta_r

  real(16),allocatable :: EX_f(:)
  complex(16),allocatable :: EX_w(:),EX_t(:)
  complex(16),allocatable :: JX_f(:),JX_w(:),JX_t(:)
  complex(16),allocatable :: GX_w(:),GX_t(:)

  integer ii,i,j,k,tmp,n

  open(101,file="./data/inverse_j.dat")
  open(102,file="./data/inverse_e.dat")
  open(103,file="./data/inverse_g.dat")
  open(111,file="./data/spectrum_j.dat")
  open(112,file="./data/spectrum_e.dat")
  open(113,file="./data/spectrum_g.dat")
  open(114,file="./data/anal2.dat")
  open(115,file="./data/anal1.dat")
  open(116,file="./data/check1.dat")
  open(117,file="./data/check2.dat")
  open(118,file="./data/check3.dat")
  open(119,file="./data/check4.dat")

  open(200,file="./data/model_env.dat")

  read(200,*) dx
  read(200,*) dy
  read(200,*) dz
  read(200,*) dtmp
  read(200,*) NN
  read(200,*) tmp
  read(200,*) tmp
  read(200,*) omega_0
  read(200,*) fmax_w
  close(200)

  allocate(EX_f(NN), EX_w(NN), EX_t(NN))
  allocate(JX_f(NN), JX_w(NN), JX_t(NN))
  allocate(GX_w(NN), GX_t(NN))

  t0      = PI/fmax_w
  beta    = PI*(fmax_w*fmax_w)
  vc      = qsqrt(2.d0*omega_0 / MU0 /1.d0)
  courant = 1.d0/vc/qsqrt(1.d0/dx/dx + 1.d0/dy/dy + 1.d0/dz/dz)
  dt      = courant

  r       = qsqrt((sx-rx)**2.d0 + (sy-ry)**2.d0 + (sz-rz)**2.d0)

  EX=0.d0
  EX_p=0.d0
  write(115,98) 0.d0,0.d0,0.d0
  do t=2,NN
    time    = t*dt
    theta   = qsqrt(MU0*sig/4.d0/time)
    theta_r = theta*r
    coef1   = moment/(4.d0*PI*sig*(r**3.d0)) * ((4.d0/qsqrt(PI)*(theta_r**3.d0)+&
    &         6.d0 /qsqrt(PI)*theta_r)*qexp(-(theta_r**2.d0)) +3.d0*erfc(theta_r))
    coef2   = moment/(4.d0*PI*sig*(r**3.d0)) * ((4.d0/qsqrt(PI)*(theta_r**3.d0)+&
    &         2.d0 /qsqrt(PI)*theta_r)*qexp(-(theta_r**2.d0)) +erfc(theta_r))
    EX      = coef1 * ((rx/r)**2.d0) - coef2*1.d0

    d_EX    = (EX-EX_p)/dt
    EX_p    = EX

    write(115,98) time,d_EX,EX
  enddo


  EX=0.d0
  do t=1,NN
    time  = t*dt
    t1    = time - r/vc
    coef1 = (MU0/4.d0/PI/r) * ((rx/r)**2.d0 -1.d0)
    coef2 = (MU0/4.d0/PI/r) * (3.d0*(rx/r)**2.d0 -1.d0)

    gamma1   = qsqrt(beta/PI)*qexp(-beta*((t1-t0)**2.d0))
    d_gamma  = -2.d0*beta*(t1-t0)*qsqrt(beta/PI)*qexp(-beta*((t1-t0)**2))
    d2_gamma = qsqrt(beta/PI)*qexp(-beta*((t1-t0)**2.d0))* &
      &        (4.d0*(beta**2.d0)*((t1-t0)**2.d0) - 2.d0*beta)

    EX = coef1*d2_gamma + coef2*(vc/r*d_gamma + (vc/r)**2.d0 * gamma1)
    write(114,*)  EX
    EX_f(t) =EX
  enddo

  open(300,file="./data/waveformX.dat")
  do i=1,NN
    read(300,*) dtmp, ReF
    JX_f(i) = ReF
  enddo

  write(*,*)  "dt      =",dt
  write(*,*)  "t0      =",t0
  write(*,*)  "omega_0 =",omega_0
  write(*,*)  "fmax_w  =",fmax_w

  do n=1,NN
    om = 2.d0*PI/NN/dt
    EX_w(n) = 0.d0
    JX_w(n) = 0.d0
    do k=1,NN
      EX_w(n) = EX_w(n) + EX_f(k) * dt* &
          &   qexp(-qsqrt(omega_0*om*(n-1))*(k-1)*dt)*cqexp(qcmplx(0.d0,1.d0)*qsqrt(omega_0*om*(n-1))*(k-1)*dt)
      JX_w(n) = JX_w(n) + JX_f(k) * dt* &
          &   qexp(-qsqrt(omega_0*om*(n-1))*(k-1)*dt)*cqexp(qcmplx(0.d0,1.d0)*qsqrt(omega_0*om*(n-1))*(k-1)*dt)
      if(n==150) write(118, 96) real(k),abs(JX_w(n)),abs(JX_f(k)),exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt),abs(JX_f(k))*exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt)
      if(n==300) write(116, 96) real(k),abs(JX_w(n)),abs(JX_f(k)),exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt),abs(JX_f(k))*exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt)
      if(n==600) write(117, 96) real(k),abs(JX_w(n)),abs(JX_f(k)),exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt),abs(JX_f(k))*exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt)
      if(n==900) write(119, 96) real(k),abs(JX_w(n)),abs(JX_f(k)),exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt),abs(JX_f(k))*exp(-sqrt(omega_0*om*(n-1))*(k-1)*dt)
    enddo
!    JX_w(n) = 2.d0*omega_0*qexp(-qsqrt(omega_0*om*(n-1))*t0) * &
!      &       cqexp(CMPLX(0.d0,1.d0)*qsqrt(omega_0*om*(n-1))*t0)* &
!      &       cqexp(-CMPLX(0.d0,1.d0)*omega_0*om*(n-1)/2.d0/beta)
!    GX_w(n) =-qcmplx(0.d0,1.d0)*qcmplx(1.d0,1.d0)*qsqrt(om*(n-1)*omega_0)* &
!      &       qexp(-qsqrt(omega_0*om*(n-1))*t0) * &
!      &       cqexp(qcmplx(0.d0,1.d0)*qsqrt(omega_0*om*(n-1))*t0)* &
!      &       cqexp(-qcmplx(0.d0,1.d0)*omega_0*om*(n-1)/2.d0/beta)
!    GX_w(n) = EX_w(n)/JX_w(n)

    write(112, 99) n*om,real(EX_w(n)),aimag(EX_w(n)),abs(EX_w(n)),atan(aimag(EX_w(n))/real(EX_w(n)))
    write(111, 99) n*om,real(JX_w(n)),aimag(JX_w(n)),abs(JX_w(n)),atan(aimag(JX_w(n))/real(JX_w(n)))
    write(113, 99) n*om,real(GX_w(n)),aimag(GX_w(n)),abs(GX_w(n)),atan(aimag(GX_w(n))/real(GX_w(n)))
  enddo

  do k=1,NN
    do n=1,NN
      JX_t(k) = JX_t(k) + JX_w(n)*cqexp(-qcmplx(0.d0,1.d0)*2.d0*PI*(k-1)*(n-1)/NN) /NN/dt
      EX_t(k) = EX_t(k) + EX_w(n)*cqexp(-qcmplx(0.d0,1.d0)*2.d0*PI*(k-1)*(n-1)/NN) /NN/dt
      GX_t(k) = GX_t(k) + GX_w(n)*cqexp(-qcmplx(0.d0,1.d0)*2.d0*PI*(k-1)*(n-1)/NN) /NN/dt
    enddo
    write(101, 98) k*dt, real(JX_t(k)), aimag(JX_t(k))
    write(102, 98) k*dt, real(EX_t(k)), aimag(EX_t(k))
    write(103, 98) k*dt, real(GX_t(k)), aimag(GX_t(k))
  enddo

99  format(E15.6, 2X, 4(E20.10,2X))
98  format(E15.6, 2X, 2(E20.10,2X))
97  format(E15.6, 4X, E50.35)
96  format(F10.3, 2X, 4(E20.10,2X))
end program quad
