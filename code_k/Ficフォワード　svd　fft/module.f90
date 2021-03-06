!-------------------------------------------------------------------------------
! モジュール
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! 物理定数
!-------------------------------------------------------------------------------
module consts
  real(8),parameter :: pi=3.141592653589793d0    ! 円周率 π
  real(8),parameter :: c=2.998d8                 ! 光速 c (m/sec)
  real(8),parameter :: eps0=8.854d-12        ! 真空の誘電率 (F/m)
  real(8),parameter :: mu0=4.0d-7*pi             ! 真空の透磁率 (H/m)
  real(8),parameter :: z0=120.0d0*pi             ! 波動インピーダンス (Ω)
end module

!-------------------------------------------------------------------------------
! FDTD module
!-------------------------------------------------------------------------------
module fdtd
  ! *** 波源情報 ***
  real(8) :: freq
  real(8) :: omg

  ! *** 時間ステップ ***
  integer :: ntime_start,ntime_end,ntime_step
  integer :: nperiod,n    ! 時間更新の間隔は(T/nperiod) ただし、Tは周期
  real(8) :: time,dt

  ! *** フィールド ***
  integer,parameter :: mx=250,my=250,mz=250
  integer :: nx,ny,nz
  real(8) :: dx,dy,dz
  integer :: nlambda0    ! 1セルの一辺の長さは(1/nlambda)
  real(8) :: ex(mx,my,mz),ey(mx,my,mz),ez(mx,my,mz)
  real(8) :: hx(mx,my,mz),hy(mz,my,mz),hz(mx,my,mz)

  ! *** 媒質定数 ***
  integer,parameter :: mmedia=10
  integer :: nmedia              ! 媒質の数
  real(8) :: eps(mmedia)         ! 誘電率 ε [F/m]
  real(8) :: mu(mmedia)          ! 透磁率 μ [H/m]
  real(8) :: sig(mmedia)         ! 導電率 σ [S/m]
  real(8) :: eps_p(mmedia)       ! ε’
  real(8) :: c_p(mmedia)         ! 仮想伝播速度



  integer :: media_id(mx+1,my+1,mz+1)

  ! *** フィールド更新係数 ***
  real(8) :: cex0,cey0,cez0
  real(8) :: cexry0,cexrz0,ceyrz0,ceyrx0,cezrx0,cezry0
  real(8) :: chxry0,chxrz0,chyrz0,chyrx0,chzrx0,chzry0
  real(8) :: cex(mmedia),cey(mmedia),cez(mmedia)
  real(8) :: cexry(mmedia),cexrz(mmedia),ceyrz(mmedia),ceyrx(mmedia),cezrx(mmedia),cezry(mmedia)
  real(8) :: cexry_p(mmedia),cexrz_p(mmedia),ceyrz_p(mmedia),ceyrx_p(mmedia),cezrx_p(mmedia),cezry_p(mmedia)
  real(8) :: chxry(mmedia),chxrz(mmedia),chyrz(mmedia),chyrx(mmedia),chzrx(mmedia),chzry(mmedia)

  ! *** Murの吸収境界条件の更新 過去の値の保存変数 ***
  real(8) :: cxd(mmedia),cyd(mmedia),czd(mmedia),cxu(mmedia),cyu(mmedia),czu(mmedia)
  real(8) :: cxx(mmedia),cyy(mmedia),czz(mmedia)
  real(8) :: cxfyd(mmedia),cxfzd(mmedia),cyfxd(mmedia),cyfzd(mmedia),czfxd(mmedia),czfyd(mmedia)
!  real(8) :: eyx1(4,my,mz),eyx2(4,my,mz),ezx1(4,my,mz),ezx2(4,my,mz)
!  real(8) :: exy1(mx,4,mz),exy2(mx,4,mz),ezy1(mz,4,mz),ezy2(mx,4,mz)
!  real(8) :: exz1(mx,my,4),exz2(mx,my,4),eyz1(mx,my,4),eyz2(mx,my,4)
!  real(8) :: hyx1(4,my,mz),hyx2(4,my,mz),hzx1(4,my,mz),hzx2(4,my,mz)
!  real(8) :: hxy1(mx,4,mz),hxy2(mx,4,mz),hzy1(mx,4,mz),hzy2(mx,4,mz)
!  real(8) :: hxz1(mz,my,4),hxz2(mx,my,4),hyz1(mx,my,4),hyz2(mx,my,4)

  real(8) :: eyx1(24,my,mz),eyx2(24,my,mz),ezx1(24,my,mz),ezx2(24,my,mz)
  real(8) :: exy1(mx,24,mz),exy2(mx,24,mz),ezy1(mz,24,mz),ezy2(mx,24,mz)
  real(8) :: exz1(mx,my,24),exz2(mx,my,24),eyz1(mx,my,24),eyz2(mx,my,24)
  real(8) :: hyx1(24,my,mz),hyx2(24,my,mz),hzx1(24,my,mz),hzx2(24,my,mz)
  real(8) :: hxy1(mx,24,mz),hxy2(mx,24,mz),hzy1(mx,24,mz),hzy2(mx,24,mz)
  real(8) :: hxz1(mz,my,24),hxz2(mx,my,24),hyz1(mx,my,24),hyz2(mx,my,24)

  ! *** PML          ***
!  integer :: i0,i1,j0,j1,k0,k1,l1,l2,l3,m,lpml
!  real(8) :: lpmlii(6,2),lpmljj(6,2),lpmlkk(6,2),R0
!  real(8) :: exy(mx,my,mz),exz(mx,my,mz),eyx(mx,my,mz),eyz(mx,my,mz),ezx(mx,my,mz),ezy(mx,my,mz) ! OK
!  real(8) :: hxy(mx,my,mz),hxz(mx,my,mz),hyx(mx,my,mz),hyz(mx,my,mz),hzx(mx,my,mz),hzy(mx,my,mz) ! OK
!  real(8) :: cxe(mx),cye(my),cze(mz)
!  real(8) :: cxh(mx),cyh(my),czh(mz)
!  real(8) :: cxel(mx),cyel(my),czel(mz)
!  real(8) :: cxhl(mx),cyhl(my),czhl(mz)
!  real(8) :: shx(mx),shy(my),shz(mz),sex(mx),sey(my),sez(mz) !OK
!  real(8) :: shmx,semx



end module
