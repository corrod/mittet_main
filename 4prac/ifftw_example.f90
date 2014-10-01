program idft

  implicit none

  integer                 :: N, N1, i, ios
  real(8)   , allocatable :: acc(:), F(:), P(:)
  complex(8), allocatable :: C(:), cacc(:)

  character(7), parameter :: inp1='fas.out', inp2='fps.out', out ='acc.out'

  integer :: plan(8)

! FFTW3を呼び出すのに必要なヘッダーファイルを include する
  include 'fftw3.f'

! ファイル（データ）の長さ N を調べる
  open(51,FILE=inp1,ACTION='read')
  N = 0
  do
    read(51,'(F12.0)',IOSTAT=ios)
    if (ios<0) exit	! ファイルの末尾に来たらループを抜ける
    N = N + 1
  end do
  close(51)

! フーリエ逆変換により作成される時刻歴データの長さ N を求め
! 加速度データ，複素フーリエ係数，フーリエ振幅，フーリエ位相の配列を確保
  N = (N-1) * 2
  allocate( acc(1:N), cacc(1:N), C(0:N-1), F(0:N/2), P(0:N/2) )

! フーリエ振幅スペクトル，フーリエ位相スペクトルを読み込む
  open(51,FILE=inp1,ACTION='read')
  read(51,'(E12.4)') F(0:N/2)
  close(51)

  open(51,FILE=inp2,ACTION='read')
  read(51,'(F12.5)') P(0:N/2)
  close(51)

! 複素フーリエ係数を求める
  C = F * CMPLX( COS(P), SIN(P) );
  C(N/2+1:N-1) = CONJG( C(N/2-1:1:-1) );

! FFTW3 を呼び出してフーリエ逆変換を行う
! 1) plan を作成する
  call dfftw_plan_dft_1d( plan, N, C, cacc, FFTW_BACKWARD, FFTW_ESTIMATE )

! 2) FFT を実行する
  call dfftw_execute(plan)

! 3) plan を破棄する
  call dfftw_destroy_plan(plan)

! フーリエ逆変換で得られたもののを規格化し，実部を取り出す
  acc = DBLE( cacc / DBLE(N) )

! 加速度時刻歴の書き出し
  open(61,FILE=out,STATUS='replace',ACTION='write')
  write(61,'(F12.5)') acc
  close(61)

end program idft
