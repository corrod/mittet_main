!!!初期値/モデル設定,変数定義!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   omega0=2πf0
!   f0=1,0Hz
!   sigmawa=3.2S/m
!   x=(i-1)*dx
!   ∂n=Σαの書き方
!   代入必要？→E'(x,omega')=E(x,omega)
!   H'(x,omega')=sqrt(-iomega/2omega0)H(x,omega)
!   J'(x,omega')=sqrt(-iomega/2omega0)J(x,omega)
!   K'(x,omega')=K(x,omega)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module const_para
    implicit none

    integer :: i,j,k
    integer, parameter :: nstep    = 900 !総タイムステップ数
    integer, parameter :: nx       = 101, ny = 101, nz = 101 !グリッド数
    integer, parameter :: x0       = 51, y0 = 51, z0 = 51  !送信源位置
    integer, parameter :: nxpml1   = 10,nypml1=10,nzpml1=10!CPMLのgrid数
    integer, parameter :: ln       = 1 !operator half rength
    real(8), parameter :: pai      = 3.14159265358979d0 !πの値
    real(8), parameter :: tau0     = 0.02d0!1.6d-4 !送信源出力時間
    real(8), parameter :: f0       = 1.0d0 !f0が小さいとdtがでかくなる
    real(8), parameter :: omega0   = 2.0d0*pai*f0 !2πf0, !ω0
    real(8), parameter :: fmax     = 25.0d0 !送信源の最大周波数
    real(8), parameter :: Glim     = 10.4d0 ! Taylor expansion参
    real(8), parameter :: dt       = 5.0d-4 !タイムステップ長 s
    real(8), parameter :: dx       = 20.0d0, dy = 20.0d0, dz = 20.0d0!dx=1.0d-2,dy=1.0d-2,dz=1.0d-2
    real(8), parameter :: CC       = 2.997924580d0 !光速
    real(8), parameter :: sigmaair = 0.0d0     !空気の導電率 S/m
    real(8), parameter :: sigmafe  = 1.03d7 !鉄の導電率 S/m
    real(8), parameter :: sigmawa  = 3.2d0  !海水の導電率 S/m
    real(8), parameter :: myu0     = 1.2566370614d-6 !真空の透磁率 H/m
    real(8), parameter :: myurair  = 1.0d0        !空気の比透磁率
    real(8), parameter :: myurfe   = 4.0d3         !鉄の比透磁率
    real(8), parameter :: myurwa   = 0.999991d0     !海水の比透磁率
    real(8), parameter :: myuair   = myurair * myu0 !空気の透磁率 H/m
    real(8), parameter :: myufe    = myurfe * myu0   !鉄の透磁率 H/m
    real(8), parameter :: myuwa    = myurwa * myu0   !鉄の透磁率 H/m
    real(8), parameter :: epsi0    = 8.854d-12   !真空の誘電率F/m
    real(8), parameter :: epsirair = 1.0006d0 !空気の比誘電率
    real(8), parameter :: epsirfe  = 1.0d0     !鉄の比誘電率
    real(8), parameter :: epsirwa  = 81.0d0      !海水の比誘電率
    real(8), parameter :: epsiair  = epsirair * epsi0 !空気の比誘電率
    real(8), parameter :: epsife   = epsirfe * epsi0   !鉄の比誘電率
    real(8), parameter :: epsiwa   = epsirwa * epsi0   !海水の比誘電率
!     real(8)            :: sigmaxx(nx,ny,nz) !diagonal sigma x
!     real(8)            :: sigmayy(nx,ny,nz) !diagonal sigma y
!     real(8)            :: sigmazz(nx,ny,nz) !diagonal sigma z

            end module const_para
