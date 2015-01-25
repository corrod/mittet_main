!!!初期値/モデル設定,変数定義!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  optimization sheme ln=2
!  from optimization scheme
!
!   omega0=2πf0
!   f0=1,0Hz
!   sigwa=3.2S/m
!   x=(i-1)*dx
!   ∂n=Σαの書き方
!   代入必要？→E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module const_para
    implicit none

    integer :: i, j, k
    integer :: sp !source position
    character(16) :: file_sp ! = sp

    integer, parameter :: nstep = 4351 ! 総タイムステップ数 　
!     integer, parameter :: nx=101, ny=101, nz=101  !グリッド数 奇数　
!     integer, parameter :: nx=101, ny=61, nz=71  !グリッド数 奇数　
    integer, parameter :: nx=131, ny=51, nz=51  !グリッド数 奇数　
    real(8), parameter :: dx = 4.5d-3, dy=dx, dz=dx !　
    real(8), parameter :: fmax = 2.0d3!1.0d2 !25.0d0 !12.5kusuda!送信源の最大周波数 　
    integer, parameter :: ncpml = 10 !CPMLのgrid数

    real(8) :: t_res(nstep) !residual time
    real(8) :: signal_res(nstep) !residual_reverse receiver

    real(8), parameter :: pi = 3.14159265358979d0 !πの値
    real(8), parameter :: f0 = 1.0d0 !f0が小さいとdtがでかくなる
    real(8), parameter :: omega0 = 2.0d0*pi*f0 !2πf0, !ω0
    complex(kind(0d0)),parameter :: I_u =(0.0d0,1.0d0)  !imaginary unit

    !optimized L=2
!     real(8), parameter :: Glim = 6.70d0
!     real(8), parameter :: c1 = 1.14443d0,c2 = - 0.04886d0
    !taylor L=2
    real(8), parameter :: Glim = 10.4d0
    real(8), parameter  :: c1 = 1.125d0, c2 = - 0.04167d0


!model
    integer, parameter :: plate = ncpml + 6 !16(15)mm   鉄板厚さ
    integer, parameter :: offset = 9        !30mm + 15mm = 45mm  送信源ープレート間の距離
    integer, parameter :: L_ver = 6         !31.1(30)mm パターン１(２)の受信間距離(縦)
    integer, parameter :: L_hori = 8        !40mm       パターン１(２)の受信間距離(横)
    integer, parameter :: L_sr = 3          !16.3(15)mm パターン２の送受信距離(縦)
    integer, parameter :: x0 = (nx+1)/2, y0 = (ny+1)/2, z0 = (nz+1)/2 !中心位置

    integer :: x_source, y_source, z_source !送信源1
    integer :: x_source2, y_source2, z_source2!送信源2
    !     integer, parameter :: x_source = x0, y_source = y0, z_source = plate + offset !送信源位置1
    !     integer, parameter :: x_source2 = x0 + L_hori, y_source2 = y0, z_source2 = plate + offset !送信源位置2

        !パターン1　ソース位置＝レシーバ①
    !     integer, parameter :: x1 = x_source,    y1 = y_source, z1 = z_source  !レシーバ位置①
    !     integer, parameter :: x2 = x1,          y2 = y1,       z2 = z1 - L_ver      !レシーバ位置②
    !     integer, parameter :: x3 = x1 + L_hori, y3 = y1,       z3 = z1 - L_ver      !レシーバ位置③
        !①
        !●
        !② ③

    !パターン2　ソース位置＝レシーバ①②間 source_posi.f90
    integer :: xx1, yy1, zz1 !レシーバ位置①
    integer :: xx2, yy2, zz2 !レシーバ位置②
    integer :: xx3, yy3, zz3 !レシーバ位置③
    !     integer, parameter :: xx1 = x_source,     yy1 = y_source, zz1 = z_source + L_sr  !レシーバ位置①
    !     integer, parameter :: xx2 = xx1,          yy2 = yy1,      zz2 = z_source - L_sr  !レシーバ位置②
    !     integer, parameter :: xx3 = xx1 + L_hori, yy3 = yy1,      zz3 = z_source - L_sr  !レシーバ位置③
    !①
    !● ●
    !② ③


!媒質パラメータ
!conductivity 導電率
    real(8), parameter :: sigair = 1.0d-11!from muhammad  !空気の導電率 S/m
    real(8), parameter :: sigfe  = 1.0d3! 　　7.5d6  !1.03d7 !鉄の導電率 S/m
    real(8), parameter :: sigwa  = 3.2d0  !海水の導電率 S/m
    real(8), parameter :: sigmin = min(sigwa,sigfe)
    real(8), parameter :: sigmax = max(sigwa,sigfe)
!permeability 透磁率
    real(8), parameter :: MU0      = 1.2566370614d-6 !真空の透磁率 H/m
    real(8), parameter :: myurair  = 1.0d0        !空気の比透磁率
    real(8), parameter :: myurfe   = 1.0d0!　　         !鉄の比透磁率
    real(8), parameter :: myurwa   = 0.999991d0     !海水の比透磁率
    real(8), parameter :: myuair   = myurair * MU0 !空気の透磁率 H/m
    real(8), parameter :: myufe    = myurfe * MU0   !鉄の透磁率 H/m
    real(8), parameter :: myuwa    = myurwa * MU0   !鉄の透磁率 H/m
!permittivity 誘電率
    real(8), parameter :: epsi0    = 8.854d-12   !真空の誘電率F/m
    real(8), parameter :: epsirair = 1.0006d0  !空気の比誘電率
    real(8), parameter :: epsirfe  = 1.0d0     !鉄の比誘電率
    real(8), parameter :: epsirwa  = 81.0d0      !海水の比誘電率
    real(8), parameter :: epsiair  = epsirair * epsi0 !空気の比誘電率
    real(8), parameter :: epsife   = epsirfe * epsi0   !鉄の比誘電率
    real(8), parameter :: epsiwa   = epsirwa * epsi0   !海水の比誘電率

    real(8) :: sig(nx,ny,nz)
    real(8) :: myu(nx,ny,nz)


!伝播速度設定 cmax, cmin
    real(8), parameter :: CC = 2.997924580d8 !光速
    real(8), parameter :: cair = sqrt(2.0d0*omega0/myuair/sigair)
    real(8), parameter :: cwa = sqrt(2.0d0*omega0/myuwa/sigwa)
    real(8), parameter :: cfe = sqrt(2.0d0*omega0/myufe/sigfe) !myufe >> myuwa 　　
    real(8), parameter :: cmin = min(cwa,cfe)
    real(8), parameter :: cmax = max(cwa,cfe)


! タイムステップ長 dt
    ! optimized dt
    ! real(8) :: dt = dx/cmax/sqrt(3.0d0)/1.19329d0
    ! taylor dt
    real(8) :: dt = 0.999d0*dx/cmax/ sqrt(3.0d0) /1.16667d0
    ! fourier dt
    ! real(8) :: dt = 0.999d0*(2.0d0*dx)/((3.0d0**0.5d0)*pi*cmax)


!CPML
    integer, parameter  :: m         = 4, ma = 1 !mの代わりにnn
    real(8), parameter  :: kappa_max = 1.0d0 !!!
    real(8), parameter  :: a_max     = 0.2d0     !!!
    real(8), parameter  :: nn        = 3.0d0 !nn should be [2,6]
    real(8), parameter  :: order     = 0.0d0 !order should be (0,3]
    real(8), parameter  :: optToMax  = 10.0d0
    real(8), parameter  :: Rcoef     = 0.01d0 !R should be [10^-2, 10^-12]
    real(8), parameter  :: epsir     = 1.0d0
    real(8) :: sig_max
    real(8) :: delta

    real(8) :: sigmax2
    real(8) :: gradmax
    real(8) :: scaler
    real(8) :: grad(nx,ny,nz)
    real(8) :: sig2(nx,ny,nz)
    real(8) :: eps2

    real(8) :: ca_x(nx,ny,nz)
    real(8) :: ca_y(nx,ny,nz)
    real(8) :: ca_z(nx,ny,nz)
    real(8) :: da_x(nx,ny,nz)
    real(8) :: da_y(nx,ny,nz)
    real(8) :: da_z(nx,ny,nz)
    real(8) :: cb_x(nx,ny,nz)
    real(8) :: cb_y(nx,ny,nz)
    real(8) :: cb_z(nx,ny,nz)
    real(8) :: db_x(nx,ny,nz)
    real(8) :: db_y(nx,ny,nz)
    real(8) :: db_z(nx,ny,nz)

    complex(kind(0d0)) :: psi_Eyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Ezx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Ezy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Exy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Exz1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Eyz1(nx,ny,nz)

    complex(kind(0d0)) :: psi_Hyx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hzx1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hzy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hxy1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hxz1(nx,ny,nz)
    complex(kind(0d0)) :: psi_Hyz1(nx,ny,nz)

    real(8) :: be_x(nx)
    real(8) :: be_y(ny)
    real(8) :: be_z(nz)
    real(8) :: bh_x(nx)
    real(8) :: bh_y(ny)
    real(8) :: bh_z(nz)
    real(8) :: ce_x(nx)
    real(8) :: ce_y(ny)
    real(8) :: ce_z(nz)
    real(8) :: ch_x(nx)
    real(8) :: ch_y(ny)
    real(8) :: ch_z(nz)

    real(8) :: esig_x(nx)
    real(8) :: esig_y(ny)
    real(8) :: esig_z(nz)
    real(8) :: msig_x(nx)
    real(8) :: msig_y(ny)
    real(8) :: msig_z(nz)
    real(8) :: ae_x(nx)
    real(8) :: ae_y(ny)
    real(8) :: ae_z(nz)
    real(8) :: am_x(nx)
    real(8) :: am_y(ny)
    real(8) :: am_z(nz)
    real(8) :: ekappa_x(nx)
    real(8) :: ekappa_y(ny)
    real(8) :: ekappa_z(nz)
    real(8) :: mkappa_x(nx)
    real(8) :: mkappa_y(ny)
    real(8) :: mkappa_z(nz)
    real(8) :: kedx(nx)
    real(8) :: kedy(ny)
    real(8) :: kedz(nz)
    real(8) :: khdx(nx)
    real(8) :: khdy(ny)
    real(8) :: khdz(nz)


!mur 変数
    real(8) :: cxd, cxu, cxx
    real(8) :: cxfyd, cxfzd
    real(8) :: cyd, cyu, cyy
    real(8) :: cyfxd, cyfzd
    real(8) :: czd, czu, czz
    real(8) :: czfxd, czfyd
    complex(kind(0d0)) :: eyx1(nx,ny,nz),eyx2(nx,ny,nz),ezx1(nx,ny,nz),ezx2(nx,ny,nz)
    complex(kind(0d0)) :: exy1(nx,ny,nz),exy2(nx,ny,nz),ezy1(nx,ny,nz),ezy2(nx,ny,nz)
    complex(kind(0d0)) :: exz1(nx,ny,nz),exz2(nx,ny,nz),eyz1(nx,ny,nz),eyz2(nx,ny,nz)

!不明な変数
    real(8) :: epsi(nx,ny,nz)
!     real(8)            :: sigxx(nx,ny,nz) !diagonal sig x
!     real(8)            :: sigyy(nx,ny,nz) !diagonal sig y
!     real(8)            :: sigzz(nx,ny,nz) !diagonal sig z

        end module const_para
