!-------------------------------------------------------------------------------
! 初期化
!-------------------------------------------------------------------------------
subroutine init
  use consts
  use fdtd
  implicit none

  integer :: id,i,j,k
  real(8) :: lambda0

  write(*,*) "**** INIT ****"
  write(1,*) "**** INIT ****"

  ! *** ファイル入力 ***
  call input

  ! *** タイムステップ、格子間隔 ***
  ! Courant の安定条件

  lambda0=c/freq

!  dx=lambda0/nlambda0
  dx=1.0d-2
  dy=dx
  dz=dx
  
  dt=(1.0d0/freq)/nperiod
!  dt=0.99*(1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c)
!  if (dt > 1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c) then
!     print *,'Sub(init): Courant stability condition must be satisfied'
     print *,'dt=',dt
!     print *,'dt<',1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2)/c
!     stop
!  end if
  write(1,*) 'dt=' ,dt

! *** 媒質定数設定 ***
  nmedia=3
! 1: 真空
  eps(1)=eps0
  mu(1)=mu0
  sig(1)=0.0d0
! 2:PEC,PMC (これは別のルーチンで処理する）
! 3:誘電体（テフロン）
  eps(3)=81.0d0*eps0
  mu(3)=mu0
  sig(3)=3.2d0

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
      cex(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))
      cey(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))
      cez(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id))))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))

      cexry(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dy
      cexrz(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dz
      ceyrz(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dz
      ceyrx(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dx
      cezrx(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dx
      cezry(id)=(dt/eps(id))/(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))/dy

      chxry(id)=(dt/mu(id))/dy
      chxrz(id)=(dt/mu(id))/dz
      chyrz(id)=(dt/mu(id))/dz
      chyrx(id)=(dt/mu(id))/dx
      chzrx(id)=(dt/mu(id))/dx
      chzry(id)=(dt/mu(id))/dy
   end do

! *** 吸収境界条件の更新定数 ***
   cxd=(c*dt-dx)/(c*dt+dx)
   cyd=(c*dt-dy)/(c*dt+dy)
   czd=(c*dt-dz)/(c*dt+dz)

   cxu=(c*dt-dx)/(c*dt+dx)
   cyu=(c*dt-dy)/(c*dt+dy)
   czu=(c*dt-dz)/(c*dt+dz)

   cxx=(2.0d0*dx)/(c*dt+dx)
   cyy=(2.0d0*dy)/(c*dt+dy)
   czz=(2.0d0*dz)/(c*dt+dz)

   cxfyd=(dx*(c*dt)**2)/(2.0d0*((dy)**2)*(c*dt+dx))
   cxfzd=(dx*(c*dt)**2)/(2.0d0*((dz)**2)*(c*dt+dx))

   cyfxd=(dy*(c*dt)**2)/(2.0d0*((dx)**2)*(c*dt+dy))
   cyfzd=(dy*(c*dt)**2)/(2.0d0*((dz)**2)*(c*dt+dy))

   czfxd=(dz*(c*dt)**2)/(2.0d0*((dx)**2)*(c*dt+dz))
   czfyd=(dz*(c*dt)**2)/(2.0d0*((dy)**2)*(c*dt+dz))

! *** PMLの吸収境界条件の更新係数 ***

   do i=1,mx
      do j=1,my
         do k=1,mz
            exy(i,j,k)=0.0d0
            exz(i,j,k)=0.0d0
            eyx(i,j,k)=0.0d0
            eyz(i,j,k)=0.0d0
            ezx(i,j,k)=0.0d0
            ezy(i,j,k)=0.0d0
            hxy(i,j,k)=0.0d0
            hxz(i,j,k)=0.0d0
            hyx(i,j,k)=0.0d0
            hyz(i,j,k)=0.0d0
            hzx(i,j,k)=0.0d0
            hzy(i,j,k)=0.0d0
         end do
      end do
   end do

   do i=1,mx
      sex(i)=0.0d0
      shx(i)=0.0d0
   end do

   do j=1,my
      sey(j)=0.0d0
      shy(j)=0.0d0
   end do

   do k=1,mz
      sez(k)=0.0d0
      shz(k)=0.0d0
   end do


   lpml=8
   m=2
   R0=1.0d-5

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

   semx=(-(m+1)*eps0*c*log(dabs(R0)))/(2.0d0*lpml*dx)
!   semx=10.0d0
   shmx=(mu0*semx)/eps0

   print *, semx,shmx 

   do i=1,lpml
      sex(i)=semx*((real(lpml)-real(i)+1.0d0)/real(lpml))**m
!      shx(i)=shmx*((real(lpml)-real(i)+0.5d0)/real(lpml))**m
      shx(i)=(mu0*sex(i))/eps0
      sex(nx-i+1)=sex(i)
      shx(nx-i+1)=shx(i)
   end do

   do i=1,nx
      cxh(i)=(2.0d0*mu0-shx(i)*dt)/(2.0d0*mu0+shx(i)*dt)
      cxhl(i)=(2.0d0*dt)/((2.0d0*mu0+shx(i)*dt)*dx)
      cxh(nx-i)=cxh(i)
      cxhl(nx-i)=cxhl(i)
      cxe(i)=(2.0d0*eps0-sex(i)*dt)/(2.0d0*eps0+sex(i)*dt)
      cxel(i)=(2.0d0*dt)/((2.0d0*eps0+sex(i)*dt)*dx)
      cxe(nx-i+1)=cxe(i)
      cxel(nx-i+1)=cxel(i)
   end do

 do j=1,lpml
      sey(j)=(float(lpml+1-j)/lpml)**m*semx
!      shy(j)=((lpml-j+0.5d0)/lpml)**m*shmx
      shy(j)=(mu0*sey(j))/eps0
      sey(ny-j+1)=sey(j)
      shy(ny-j+1)=shy(j)
 end do

 do j=1,ny
      cyh(j)=(2.0d0*mu0-shy(j)*dt)/(2.0d0*mu0+shy(j)*dt)
      cyhl(j)=(2.0d0*dt)/((2.0d0*mu0+shy(j)*dt)*dy)
      cyh(ny-j)=cyh(j)
      cyhl(nx-j)=cyhl(j)
!      cye(j)=(2.0d0*eps0-sey(j)*dt)/(2.0d0*eps0+sey(j)*dt)
      cye(j)=(1.0d0-((sey(j)*dt)/(2.0d0*eps0)))/(1.0d0+((sey(j)*dt)/(2.0d0*eps0)))
!      cyel(j)=(2.0d0*dt)/((2.0d0*eps0+sey(j)*dt)*dy)
!cezyry(j)=(dt/eps0)/(1.0d0+((esigy(j)*dt)/(2.0d0*eps0)))/dyp

      cyel(j)=(dt/eps0)/((1.0d0+((sey(j)*dt)/(2.0d0*eps0)))*dy)
      cye(ny-j+1)=cye(j)
      cyel(ny-j+1)=cyel(j)
   end do



 do k=1,lpml
      sez(k)=(float(lpml+1-k)/lpml)**m*semx
!      shz(k)=((lpml-k+0.5d0)/lpml)**m*shmx
      shz(k)=(mu0*sez(k))/eps0
      sez(nz-k+1)=sez(k)
      shz(nz-k+1)=shz(k)
  end do

  do k=1,nz
      czh(k)=(2.0d0*mu0-shz(k)*dt)/(2.0d0*mu0+shz(k)*dt)
      czhl(k)=(2.0d0*dt)/((2.0d0*mu0+shz(k)*dt)*dz)
      czh(nz-k)=czh(k)
      czhl(nz-k)=czhl(k)
      cze(k)=(2.0d0*eps0-sez(k)*dt)/(2.0d0*eps0+sez(k)*dt)
      czel(k)=(2.0d0*dt)/((2.0d0*eps0+sez(k)*dt)*dz)
      cze(nz-k+1)=cze(k)
      czel(nz-k+1)=czel(k)
   end do
   
   open(3,file='pml_conductivity_y')
   open(31,file='ce_y')
   open(32,file='ch_y')
   do j=1,ny
      write(3,*) j,sey(j),shy(j)
      write(31,*) j,cye(j),cyel(j)
      write(32,*) j,cyh(j),cyhl(j)
   end do
   close(3)
   close(31)
   close(32)

   open(4,file='pml_conductivity_x')
   open(33,file='ce_x')
   open(34,file='ch_x')
   do i=1,nx
      write(4,*) i,sex(i),shx(i)
      write(33,*) i,cye(i),cyel(i)
      write(34,*) i,cyh(i),cyhl(i)
   end do
   close(4)
   close(33)
   close(34)

   open(5,file='pml_conductivity_z')
   open(35,file='ce_z')
   open(36,file='ch_z')
   do k=1,nz
      write(5,*) k,sez(k),shz(k)
      write(35,*) k,cze(k),czel(k)
      write(36,*) k,czh(k),czhl(k)
   end do
   close(5)
   close(35)
   close(36)

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
  read(fp,*) text,nx
  read(fp,*) text,nlambda0
  read(fp,*) text,nperiod
  read(fp,*) text,ntime_start
  read(fp,*) text,ntime_end
  read(fp,*) text,ntime_step
  close(fp)

!  freq=freq*1.0d9          ! GHz → HZに変換
  ny=nx
  nz=nx
  ! エラー処理
  if((nx>mx).or.(ny>my).or.(nz>mz)) then
     write(*,*) " Sub(input): mx,my,mzare small!"
     stop
  end if

  ! パラメータ確認
  write(1,*) "FREQUENCY (GHZ)=",freq
  write(1,*) "NUMBER OF CELLS [NX,NY,NZ]=",nx
  write(1,*) "CELL INTERVAL [(FREE WAVELENGTH)/N]=",nlambda0
  write(1,*) "TIME STEP [START] (sec)=",ntime_start
  write(1,*) "TIME STEP [END]   (sec)=",ntime_end
  write(1,*) "TIME STEP [STEP]  (sec)=",ntime_step

  return
end subroutine
