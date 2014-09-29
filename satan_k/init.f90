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
  nmedia=3
! 1: 真空
  eps(1)=epsilon0
  mu(1)=mu0
  sig(1)=0.0d0
! 2:PEC,PMC (これは別のルーチンで処理する）
! 3:誘電体（海水）
  eps(3)=81.0d0*epsilon0
  mu(3)=mu0
  sig(3)=3.2d0

! *** fictitious domain ***
  do id=3,nmedia
  omg=2.0d0*pi*freq
  eps_p(id)=sig(id)/(2.0d0*omg)
  c_p(id)=sqrt((2.0d0*omg)/(mu0*sig(id)))

  print *, "omg=",omg
  print *, "eps_p(id)=",eps_p(id)
  print *, "c_p(id)=", c_p(id)

! *** タイムステップ、格子間隔 ***
! Courant の安定条件

  lambda0=c_p(id)/freq
  end do

!  dx=lambda0/nlambda0
  dx=1.0d2
  dy=dx
  dz=dx
!  dt=(1.0d0/freq)/nperiod
  dt=0.99*(1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(3))
  if (dt > 1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c_p(3)) then
     print *,'Sub(init): Courant stability condition must be satisfied'
     print *,'dt=',dt
     print *,'dt<',1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2)/c_p(3)
     stop
  end if
  write(1,*) 'dt=' ,dt


! *** フィールド更新係数 ***
! 1: 真空
  cex0=1.0d0
  cey0=1.0d0
  cez0=1.0d0

  cexry0=(dt/epsilon0)/dy
  cexrz0=(dt/epsilon0)/dz
  ceyrz0=(dt/epsilon0)/dz
  ceyrx0=(dt/epsilon0)/dx
  cezrx0=(dt/epsilon0)/dx
  cezry0=(dt/epsilon0)/dy

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
   cxd=(c_p(3)*dt-dx)/(c_p(3)*dt+dx)
   cyd=(c_p(3)*dt-dy)/(c_p(3)*dt+dy)
   czd=(c_p(3)*dt-dz)/(c_p(3)*dt+dz)

   cxu=(c_p(3)*dt-dx)/(c_p(3)*dt+dx)
   cyu=(c_p(3)*dt-dy)/(c_p(3)*dt+dy)
   czu=(c_p(3)*dt-dz)/(c_p(3)*dt+dz)

   cxx=(2.0d0*dx)/(c_p(3)*dt+dx)
   cyy=(2.0d0*dy)/(c_p(3)*dt+dy)
   czz=(2.0d0*dz)/(c_p(3)*dt+dz)

   cxfyd=(dx*(c_p(3)*dt)**2)/(2.0d0*((dy)**2)*(c_p(3)*dt+dx))
   cxfzd=(dx*(c_p(3)*dt)**2)/(2.0d0*((dz)**2)*(c_p(3)*dt+dx))

   cyfxd=(dy*(c_p(3)*dt)**2)/(2.0d0*((dx)**2)*(c_p(3)*dt+dy))
   cyfzd=(dy*(c_p(3)*dt)**2)/(2.0d0*((dz)**2)*(c_p(3)*dt+dy))

   czfxd=(dz*(c_p(3)*dt)**2)/(2.0d0*((dx)**2)*(c_p(3)*dt+dz))
   czfyd=(dz*(c_p(3)*dt)**2)/(2.0d0*((dy)**2)*(c_p(3)*dt+dz))

! *** PML吸収境界条件の更新係数 ***

   do i=1,nx
      cxh(i)=1.0d0
      cxhl(i)=dt/mu0/dx
      cxe(i)=1.0d0
      cxel(i)=dt/eps_p(3)/dx
   end do

   do j=1,ny
      cyh(j)=1.0d0
      cyhl(j)=dt/mu0/dy
      cye(j)=1.0d0
      cyel(j)=dt/eps_p(3)/dy
   end do

   do k=1,nz
      czh(k)=1.0d0
      czhl(k)=dt/mu0/dz
      cze(k)=1.0d0
      czel(k)=dt/eps_p(3)/dz
   end do


   lpml=16
   m=4
   R0=1.0d-6

   lpmlii(1,1)=1
   lpmljj(1,1)=1
   lpmlkk(1,1)=1
   lpmlii(1,2)=lpml+1
   lpmljj(1,2)=ny
   lpmlkk(1,2)=nz
   lpmlii(2,1)=nx-lpml
   lpmljj(2,1)=1
   lpmlkk(2,1)=1
   lpmlii(2,2)=nx
   lpmljj(2,2)=ny
   lpmlkk(2,2)=nz
   lpmlii(3,1)=1
   lpmljj(3,1)=1
   lpmlkk(3,1)=1
   lpmlii(3,2)=nx
   lpmljj(3,2)=lpml+1
   lpmlkk(3,2)=nz
   lpmlii(4,1)=1
   lpmljj(4,1)=ny-lpml
   lpmlkk(4,1)=1
   lpmlii(4,2)=nx
   lpmljj(4,2)=ny
   lpmlkk(4,2)=nz
   lpmlii(5,1)=1
   lpmljj(5,1)=1
   lpmlkk(5,1)=1
   lpmlii(5,2)=nx
   lpmljj(5,2)=ny
   lpmlkk(5,2)=lpml+1
   lpmlii(6,1)=1
   lpmljj(6,1)=1
   lpmlkk(6,1)=nz-lpml
   lpmlii(6,2)=nx
   lpmljj(6,2)=ny
   lpmlkk(6,2)=nz


   lpmlst(1)=1
   lpmlst(2)=lpmlst(1)+ny*nz*(lpml+1)
   lpmlst(3)=1
   lpmlst(4)=lpmlst(3)+nx*nz*(lpml+1)
   lpmlst(5)=1
   lpmlst(6)=lpmlst(5)+nx+ny*(lpml+1)

   semx=(-(m+1)*eps_p(3)*c_p(3)*log(dabs(R0)))/(2.0d0*lpml*dx)
   shmx=(mu0*semx)/eps(3)

   print *, semx,shmx

   do i=1,lpml
      sex=(float(lpml+1-i)/lpml)**m*semx
      shx=((lpml-i+0.5d0)/lpml)**m*shmx
      cxh(i)=(2.0d0*mu0-shx*dt)/(2.0d0*mu0+shx*dt)
      cxhl(i)=(2.0d0*dt)/((2.0d0*mu0+shx*dt)*dx)
      cxh(nx-i)=cxh(i)
      cxhl(nx-i)=cxhl(i)
      cxe(i)=(2.0d0*eps(3)-sex*dt)/(2.0d0*eps(3)+sex*dt)
      cxel(i)=(2.0d0*dt)/((2.0d0*eps(3)+sex*dt)*dx)
      cxe(nx-i+1)=cxe(i)
      cxel(nx-i+1)=cxel(i)
   end do

 do j=1,lpml
      sey=(float(lpml+1-j)/lpml)**m*semx
      shy=((lpml-i+0.5d0)/lpml)**m*shmx
      cyh(j)=(2.0d0*mu0-shy*dt)/(2.0d0*mu0+shy*dt)
      cyhl(j)=(2.0d0*dt)/((2.0d0*mu0+shy*dt)*dy)
      cyh(ny-j)=cyh(j)
      cyhl(nx-j)=cyhl(j)
      cye(j)=(2.0d0*eps(3)-sey*dt)/(2.0d0*eps(3)+sey*dt)
      cyel(j)=(2.0d0*dt)/((2.0d0*eps(3)+sey*dt)*dy)
      cye(ny-j+1)=cye(j)
      cyel(nx-j+1)=cyel(j)
   end do

 do k=1,lpml
      sez=(float(lpml+1-k)/lpml)**m*semx
      shz=((lpml-i+0.5d0)/lpml)**m*shmx
      czh(k)=(2.0d0*mu0-shz*dt)/(2.0d0*mu0+shz*dt)
      czhl(k)=(2.0d0*dt)/((2.0d0*mu0+shz*dt)*dz)
      czh(nz-k)=czh(k)
      czhl(nz-k)=czhl(k)
      cze(k)=(2.0d0*eps(3)-sez*dt)/(2.0d0*eps(3)+sez*dt)
      czel(k)=(2.0d0*dt)/((2.0d0*eps(3)+sez*dt)*dz)
      cze(nx-k+1)=cze(k)
      czel(nx-k+1)=czel(k)
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

  open(fp,file='input.dat')
  read(fp,*) text
  read(fp,*) text,freq
  read(fp,*) text,nx,ny,nz
  read(fp,*) text,nlambda0
  read(fp,*) text,nperiod
  read(fp,*) text,ntime_start
  read(fp,*) text,ntime_end
  read(fp,*) text,ntime_step
  close(fp)

!  freq=freq*1.0d9          ! GHz → HZに変換
!  ny=nx
!  nz=nx
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
