!****************************************************************************
! モジュール
!****************************************************************************

!============================================================================
! 物理定数
!============================================================================
module consts
   real(8),parameter :: pi=3.141592653589793d0        ! 円周率 π
   real(8),parameter :: c=2.998d8                     ! 光速 c (m/sec)
   real(8),parameter :: epsilon0=8.854d-12            ! 真空の誘電率 (F/m)
   real(8),parameter :: mu0=4.0d-7*pi                 ! 真空の透磁率 (H/m)
   real(8),parameter :: z0=120.0d0*pi                 ! 波動インピーダンス (Ω)
end module

!============================================================================
! FDTD module
!============================================================================
module fdtd
   ! ******** 波源情報 ********
   real(8) :: freq

   ! ******** 時間ステップ ********
   integer :: ntime_start,ntime_end,ntime_step
   integer :: nperiod    ! 時間更新の間隔は(T/nperiod) ただし、Tは周期
   real(8) :: time,dt

   ! ******** フィールド ********
   integer,parameter :: mx=100,my=100,mz=100
   integer :: nx,ny,nz
   real(8) :: dx,dy,dz
   integer :: nlambda0    ! １セルの一辺の長さは(1/nlambda)
   real(8) :: ex(mx,my,mz),ey(mx,my,mz),ez(mx,my,mz)
   real(8) :: hx(mx,my,mz),hy(mx,my,mz),hz(mx,my,mz)

   ! ******** 媒質定数 ********
   integer,parameter :: mmedia=10
   integer :: nmedia         ! 媒質の数
   real(8) :: eps(mmedia)    ! 誘電率ε [F/m]
   real(8) :: mu(mmedia)     ! 透磁率μ [H/m]
   real(8) :: sig(mmedia)    ! 導電率σ [S/m]

   integer :: media_id(mx+1,my+1,mz+1)

   ! ******** フィールド更新係数 ********
   real(8) :: cex0,cey0,cez0
   real(8) :: cexry0,cexrz0, &
              ceyrz0,ceyrx0, &
              cezrx0,cezry0
   real(8) :: chxry0,chxrz0, &
              chyrz0,chyrx0, &
              chzrx0,chzry0
   real(8) :: cex(mmedia),cey(mmedia),cez(mmedia)
   real(8) :: cexry(mmedia),cexrz(mmedia), &
              ceyrz(mmedia),ceyrx(mmedia), &
              cezrx(mmedia),cezry(mmedia)
   real(8) :: chxry(mmedia),chxrz(mmedia), &
              chyrz(mmedia),chyrx(mmedia), &
              chzrx(mmedia),chzry(mmedia)

   ! ******** 吸収境界条件の更新・過去の値の保存変数 ********
   real(8) :: cxd,cyd,czd, &
              cxu,cyu,czu
   real(8) :: cxx,cyy,czz
   real(8) :: cxfyd,cxfzd, &
              cyfxd,cyfzd, &
              czfxd,czfyd
   real(8) :: eyx1(4,my,mz),eyx2(4,my,mz),ezx1(4,my,mz),ezx2(4,my,mz), &
              exy1(mx,4,mz),exy2(mx,4,mz),ezy1(mx,4,mz),ezy2(mx,4,mz), &
              exz1(mx,my,4),exz2(mx,my,4),eyz1(mx,my,4),eyz2(mx,my,4)
   real(8) :: hyx1(4,my,mz),hyx2(4,my,mz),hzx1(4,my,mz),hzx2(4,my,mz), &
              hxy1(mx,4,mz),hxy2(mx,4,mz),hzy1(mx,4,mz),hzy2(mx,4,mz), &
              hxz1(mx,my,4),hxz2(mx,my,4),hyz1(mx,my,4),hyz2(mx,my,4)
end module

!
! End of file
!

