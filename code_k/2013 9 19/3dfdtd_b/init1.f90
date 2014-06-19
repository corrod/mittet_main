!****************************************************************************
! 初期化
!****************************************************************************
subroutine init
   use consts
   use fdtd
   implicit none

   integer :: id
   real(8) :: lambda0

   write(*,*) "**** INIT ****"
   write(1,*) "**** INIT ****"

   ! ******** ファイル入力 ********
   call input

   ! ******** タイムステップ、格子間隔 ********
   ! Courant の安定条件
   ! v dt < 1/sqrt((1/dx)**2+(1/dy)**2+(1/dz)**2)
   ! v dt < dx/sqrt(3) when dx=dy=dz
   lambda0=c/freq

   dx=lambda0/nlambda0
   dy=dx
   dz=dx
   dt=(1.0d0/freq)/nperiod
   if(dt > 1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2+(1.0d0/dz)**2)/c) then
      print *,'Sub(init): Courant stability condition must be satisfied'
      print *,'dt = ',dt
      print *,'dt < ',1.0d0/sqrt((1.0d0/dx)**2+(1.0d0/dy)**2)/c
      stop
   end if
   write(1,*) 'dt=',dt

   ! ******** 媒質定数設定 ********
   nmedia=3
   ! 1: 真空
   eps(1)=epsilon0
   mu(1)=mu0
   sig(1)=0.0d0
   ! 2: PEC,PMC (これは別のルーチンで処理する)
   ! 3: 誘電体（テフロン）
   eps(3)=2.17d0*epsilon0
   mu(3)=mu0
   sig(3)=0.0d0

   ! ******** フィールド更新係数 ********
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
      cex(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id)))) &
              /(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))
      cey(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id)))) &
              /(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))
      cez(id)=(1.0d0-((sig(id)*dt)/(2.0d0*eps(id)))) &
              /(1.0d0+((sig(id)*dt)/(2.0d0*eps(id))))

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

   ! ******** 吸収境界条件の更新定数 ********
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

   return
end subroutine

!----------------------------------------------------------------------------
! ファイル入力
!----------------------------------------------------------------------------
subroutine input
   use fdtd
   implicit none

   integer :: fp    ! ファイル識別子
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

   freq=freq*1.0d9    ! GHz → Hz に変換
   ny=nx
   nz=nx
   ! エラー処理
   if((nx>mx).or.(ny>my).or.(nz>mz)) then
      write(*,*) "Sub(input): mx, my or mz are small!"
      stop
   end if

   ! パラメータ確認
   write(1,*) "FREQUENCY (GHz)=",freq
   write(1,*) "NUMBER OF CELLS [NX,NY,NZ]=",nx
   write(1,*) "CELL INTERVAL [(FREE WAVELENGTH)/N]=",nlambda0
   write(1,*) "TIME INTERVAL [(PERIOD)/N]=",nperiod
   write(1,*) "TIME STEP [START] (sec)=",ntime_start
   write(1,*) "TIME STEP [END]   (sec)=",ntime_end
   write(1,*) "TIME STEP [STEP]  (sec)=",ntime_step

   return
end subroutine

!
! End of file
!
