!-------------------------------------------------------------------------------
! 初期化
!-------------------------------------------------------------------------------
subroutine init
  use consts
  use fdtd
  implicit none

  integer :: id,i,j,k,l
  real(8) :: lambda0

  write(*,*) "**** INIT ****"
  write(1,*) "**** INIT ****"
! *** ファイル入力 ***
  call input

! *** 電場と磁場の初期化
  do i=1,nx
     do j=1,ny
        do k=1,nz
           ex(i,j,k)=0.0d0
           ey(i,j,k)=0.0d0
           ez(i,j,k)=0.0d0
           hx(i,j,k)=0.0d0
           hy(i,j,k)=0.0d0
           hz(i,j,k)=0.0d0
        end do
     end do
  end do

! *** 媒質定数設定 ***
  nmedia=7
! 1: 真空
  eps(1)=eps0
  mu(1)=mu0
  sig(1)=0.0d0
! 2:PEC,PMC (これは別のルーチンで処理する）
! 3:誘電体（海水）
  eps(3)=81.0d0*eps0
  mu(3)=mu0
  sig(3)=3.2d0

! 4:ターゲット（高非抵抗帯）
  eps(4)=58.0d0*eps0
  mu(4)=mu0
  sig(4)=2.0d-2

! 5:誘電体（堆積層1)
  eps(5)=30.0d0*eps0
  mu(5)=mu0
  sig(5)=1.0d0

! 6:誘電体(堆積層2)
  eps(6)=30.0d0*eps0
  mu(6)=mu0
  sig(6)=0.5d0

! 7:メタンハイドレート
  eps(7)=58.0d0*eps0
  mu(7)=mu0
!  sig(7)=1.0d0/3.0d1
  sig(7)=1.0d-1

! *** fictitious domain ***
  do id=3,nmedia
  omg=2.0d0*pi*freq
  eps_p(id)=sig(id)/(2.0d0*omg)
  c_p(id)=sqrt((2.0d0*omg)/(mu0*sig(id)))


! *** タイムステップ、格子間隔 ***
! Courant の安定条件

  lambda0=c_p(id)/freq
  end do


  print *, "omg=",omg
  print *, "eps_p(id)=",eps_p(3),eps_p(4),eps_p(5),eps_p(6),eps_p(7)
  print *, "c_p(id)=", c_p(3),c_p(4),c_p(5),c_p(6),c_p(7)

!  dx=lambda0/nlambda0
  dx=2.0d1
  dy=2.0d1
  dz=2.0d1
!  dt=(1.0d0/freq)/nperiod
  dt=0.99*(1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(7))

  if (dt > 1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(7)) then
     print *,'Sub(init): Courant stability condition must be satisfied'
     print *,'dt=',dt
     print *,'dt<',1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2)/c_p(7)
     stop
  end if
  write(1,*) 'dt=' ,dt


! *** フィールド更新係数 ***
! 1: 真空
  cex0=1.0d0
  cey0=1.0d0
  cez0=1.0d0

  cexry0=(dt/eps0)/dy
  cexrz0=(dt/eps0)/dz
  ceyrz0=(dt/eps0)/dz
  ceyrx0=(dt/eps0)/dx
  cezrx0=(dt/eps0)/dx
  cezry0=(dt/eps0)/dy

  chxry0=(dt/mu0)/dy
  chxrz0=(dt/mu0)/dz
  chyrz0=(dt/mu0)/dz
  chyrx0=(dt/mu0)/dx
  chzrx0=(dt/mu0)/dx
  chzry0=(dt/mu0)/dy

! 3以上: 一般の媒質
   do id=3,nmedia
      cex(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps_p(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))
      cey(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps_p(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))
      cez(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps_p(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))

      cexry(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dy
      cexrz(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dz
      ceyrz(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dz
      ceyrx(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dx
      cezrx(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dx
      cezry(id)=(dt/eps_p(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps_p(id))))/dy

      cexry_p(id)=(dt/eps_p(id))/dy
      cexrz_p(id)=(dt/eps_p(id))/dz
      ceyrz_p(id)=(dt/eps_p(id))/dz
      ceyrx_p(id)=(dt/eps_p(id))/dx
      cezrx_p(id)=(dt/eps_p(id))/dx
      cezry_p(id)=(dt/eps_p(id))/dy

      chxry(id)=(dt/mu(id))/dy
      chxrz(id)=(dt/mu(id))/dz
      chyrz(id)=(dt/mu(id))/dz
      chyrx(id)=(dt/mu(id))/dx
      chzrx(id)=(dt/mu(id))/dx
      chzry(id)=(dt/mu(id))/dy
   end do

! *** Murの吸収境界条件の更新定数 ***
   do id=3,nmedia
   cxd(id)=(c_p(id)*dt-dx)/(c_p(id)*dt+dx)
   cyd(id)=(c_p(id)*dt-dy)/(c_p(id)*dt+dy)
   czd(id)=(c_p(id)*dt-dz)/(c_p(id)*dt+dz)

   cxu(id)=(c_p(id)*dt-dx)/(c_p(id)*dt+dx)
   cyu(id)=(c_p(id)*dt-dy)/(c_p(id)*dt+dy)
   czu(id)=(c_p(id)*dt-dz)/(c_p(id)*dt+dz)

   cxx(id)=(2.0d0*dx)/(c_p(id)*dt+dx)
   cyy(id)=(2.0d0*dy)/(c_p(id)*dt+dy)
   czz(id)=(2.0d0*dz)/(c_p(id)*dt+dz)

   cxfyd(id)=(dx*(c_p(id)*dt)**2)/(2.0d0*((dy)**2)*(c_p(id)*dt+dx))
   cxfzd(id)=(dx*(c_p(id)*dt)**2)/(2.0d0*((dz)**2)*(c_p(id)*dt+dx))

   cyfxd(id)=(dy*(c_p(id)*dt)**2)/(2.0d0*((dx)**2)*(c_p(id)*dt+dy))
   cyfzd(id)=(dy*(c_p(id)*dt)**2)/(2.0d0*((dz)**2)*(c_p(id)*dt+dy))

   czfxd(id)=(dz*(c_p(id)*dt)**2)/(2.0d0*((dx)**2)*(c_p(id)*dt+dz))
   czfyd(id)=(dz*(c_p(id)*dt)**2)/(2.0d0*((dy)**2)*(c_p(id)*dt+dz))
   end do

   return
end subroutine

!------------------------------------------------------------------
! ファイル入力
!------------------------------------------------------------------
subroutine input
  use fdtd
  implicit none

  integer :: fp            ! ファイル識別子
  character(128) :: text

  open(fp,file='input_target.dat')
  read(fp,*) text
  read(fp,*) text,freq
  read(fp,*) text,nx,ny,nz
  read(fp,*) text,nlambda0
  read(fp,*) text,nperiod
  read(fp,*) text,ntime_start
  read(fp,*) text,ntime_end
  read(fp,*) text,ntime_step
  close(fp)

  ! エラー処理
  if((nx>mx).or.(ny>my).or.(nz>mz)) then
     write(*,*) " Sub(input): mx,my,mzare small!"
     stop
  end if

  ! パラメータ確認
  write(1,*) "FREQUENCY (GHZ)=",freq
  write(1,*) "NUMBER OF CELLS [NX,NY,NZ]=",nx,ny,nz
  write(1,*) "CELL INTERVAL [(FREE WAVELENGTH)/N]=",nlambda0
  write(1,*) "TIME STEP [START] (sec)=",ntime_start
  write(1,*) "TIME STEP [END]   (sec)=",ntime_end
  write(1,*) "TIME STEP [STEP]  (sec)=",ntime_step

  return
end subroutine
