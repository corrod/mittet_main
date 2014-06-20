!!!dt,dx,dy,dzの設定cmax,cminの計算************************************************
! subroutine set_d_txyz(cmax)
subroutine set_d_txyz
    use const_para
    implicit none

    integer :: Nt1,Nt2
    real(8) :: S !タイムステップ数を求める際の係数
!     real(8) :: cwa, cfe, cair
    real(8) :: courant
!     real(8), intent(out) :: cmax
!     real(8) :: cmin
    real(8) :: t_cal !計測時間
    real(8) :: dt_max
    real(8) :: fmax_w !最大周波数（上限）
    real(8) :: dt_ideal !クーラン条件を満たすdt
    real(8) :: cfl_limit
    real(8) :: fourier_limit
    real(8) :: dt_mittet



!媒質中の伝播速度cmax,cmin計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cwa  = sqrt(2.0d0*omega0/myuwa/sigwa)   !!!!cmin=sqrt(2.0d0*omega0/myuwa/maxv)
!     cfe  = sqrt(2.0d0*omega0/myufe/sigfe)   !!!!cmax=sqrt(2.0d0*omega0/myufe/maxv)
!     cair = 1.0d0 / sqrt(myuair*epsiair)

!     cmax = cwa  !max(cwa,cfe,cair)!   最大伝播速度cmax計算
!     cmin = cwa  !min(cwa,cfe,cair)
!     cmin = sqrt(2.0d0*omega0/MU0/1.0d0)  !!別の求め方もあるだり



!計算条件の出力!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*)  'nstep',nstep
    write(*,*)  'dt',dt
    write(*,*)  'fmax', fmax
    write(*,*)  'omega0',omega0
    write(*,'(a,3i5)')      'nx,ny,nz',nx,ny,nz
    write(*,'(a,3e12.4)')   'dx,dy,dz',dx,dy,dz
    write(*,*)              '波長λ', cmax/fmax
    write(*,*) 'cwa',cwa
    write(*,*) 'cfe',cfe
    write(*,*) 'cmax',cmax
!     write(*,*) 'cmin',cmin
    write(*,*) 'propagate distance c*t', cwa*dt*nstep
    write(*,*) '１/2辺の長さdx*nx/2',dx*nx*0.5d0
    write(*,*) 'dt*nstep',dt*nstep
    write(*,*) 'sig_wa,myu_wa',sigwa,myuwa
    write(*,*) 'sig_fe,myu_fe',sigfe,myufe


!!!グリッド間隔dxの設定;グリッド分散!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(dx>(cmax/fmax)) then
    write(*,*) '*****grid dispersion may happen*****'
endif

!!!タイムステップdtの計算;クーラン条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !クーラン条件
    courant = 1.0d0/cmax/sqrt(1.0d0/dx**2.0d0 + 1.0d0/dy**2.0d0 + 1.0d0/dz**2.0d0)
    dt_max = (2.0d0*dx)/((3.0d0**0.5d0)*pi*cmax) !こっちのほうがぽい

    dt_ideal = courant*6.0d0/7.0d0*0.999d0        !デカすぎ

    dt_mittet = dx/(3*0.5d0)/cmax !standard 2nd order space schme

    write(*,*) '#courant imamu dt_max', dt_max
    write(*,*) '#courant dt_ideal',dt_ideal  !デカすぎ
    write(*,*) '#dt mittet',dt_mittet


!CFL limit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cfl_limit = cmax*dt/dx*(3.0d0**0.5d0)
    if(cfl_limit>1) then
        write(*,*) '*****CFL limit is violated*****'
    endif


!Fourier limit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fourier_limit = cmax*dt/dx*pi/2.0d0*(3.0d0**0.5d0)
    if(fourier_limit>1) then
        write(*,*) '*****Fourier limit is violated*****'
    endif


!!!the number of timestep nstepの計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     S=0.6d0  !S =[0.5:1.0]

     Nt1=int(S*nx/sqrt( 1.0d0/dx**2 +1.0d0/dy**2 + 1.0d0/dz**2 ) )

     Nt2=int(S*nx/(3.0d0**0.5))!*sqrt(sigmax/sigmin))

    write(*,*) '#number of timestep Nt1,Nt2',Nt1,Nt2

     if(Nt1 >nstep) then
    write(*,*) "******* time step nstep is violated *******"
    endif


!最大の周波数の確認!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cmin = sqrt(2.0d0*omega0/MU0/1.0d0)
    fmax_w = cmin /Glim /max(dx,dy,dz)

    write(*,*) "#fmax_w",fmax_w
     if(fmax_w <fmax) then
    write(*,*), "******* fmax is violated *******\n"
    endif
            end subroutine set_d_txyz