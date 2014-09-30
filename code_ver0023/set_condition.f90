
!/////////////////////////////////////////////////////////////////////////////////
!!!dt,dx,dy,dzの設定cmax,cminの計算
!////////////////////////////////////////////////////////////////////////////////
subroutine confirm_parameter
    use const_para
    implicit none

    integer :: Nt1,Nt2,Nt3,Nt4
    real(8) :: S !タイムステップ数を求める際の係数
!     real(8) :: cwa, cfe, cair
    real(8) :: courant
!     real(8), intent(out) :: cmax
!     real(8) :: cmin
    real(8) :: t_cal !計測時間
    real(8) :: dt2,dt4,dt_wh
    real(8) :: fmax_w !最大周波数（上限）
    real(8) :: dt_ideal !クーラン条件を満たすdt
    real(8) :: cfl_limit
    real(8) :: fourier_limit
    real(8) :: skinwa, skinfe



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
    write(*,*)              '波長λ wa, fe', cmax/fmax, cmin/fmax
!     write(*,*) 'cwa',cwa
!     write(*,*) 'cfe',cfe
    write(*,*) 'cmax=cwa,cmin=cfe',cmax,cmin
!     write(*,*) 'c_test myufe>myuwa',sqrt(2.0d0*omega0/myuair/sigfe)
    write(*,*) '反射波の到達時間 (nx-10)*dx/c', (nx-10)*dx/cmax
    write(*,*) '到達距離 cmax*t,cmin*t', cmax*dt*nstep, cmin*dt*nstep
    write(*,*) '１/2辺の長さdx*nx/2',dx*nx*0.5d0
    write(*,*) 'dt*nstep',dt*nstep
    write(*,*) 'sigfe, sigwa', sigfe, sigwa
    write(*,*) 'myufe, myuwa', myufe, myuwa
    write(*,*) 'epsife, epsiwa', epsife, epsiwa



!!!skin depth の計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    skinwa = sqrt(1.0d0 / (pi*fmax*sigwa*myuwa))
    skinfe = sqrt(1.0d0 / (pi*fmax*sigfe*myufe))

    write(*,*) '# skindepth wa, fe', skinwa, skinfe



!!! dx グリッド間隔:dxの設定;グリッド分散!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if(dx>(cmax/fmax)) then
!     write(*,*) '*****grid dispersion may happen** dx>(cmax/fmax)***************'
!     write(*,*) "# cmax/fmax" , cmax/fmax
! endif

! if(dx>(cmin/fmax)) then
!     write(*,*) '*****grid dispersion may happen** dx>(cmin/fmax)***************'
!     write(*,*) "# cmin/fmax" , cmin/fmax
! endif

if(dx>(0.999d0*cmin/fmax/Glim)) then
    write(*,*) '*****grid dispersion may happen***** dx>(0.999d0*cmin/fmax/Glim) ********'
endif
    write(*,*) '# dx should  less than (0.999d0*cmin/fmax/Glim)',(0.999d0*cmin/fmax/Glim)
!     write(*,*) '# dx should  less than (0.999d0*cmax/fmax/Glim)',(0.999d0*cmax/fmax/Glim)


!!! dt タイムステップdtの計算;クーラン条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !クーラン条件
    !standard 2nd order space schme
    dt2 = dx/(3.0d0*0.5d0)/cmax
    !fourth order scheme with taylor operator
    dt4 = (2.0d0*dx)/((3.0d0**0.5d0)*pi*cmax) !こっちのほうがぽい
    !参照不明
    courant = 1.0d0/cmax/sqrt(1.0d0/dx**2.0d0 + 1.0d0/dy**2.0d0 + 1.0d0/dz**2.0d0)
    dt_ideal = courant*0.999d0!6.0d0/7.0d0*0.999d0        !デカすぎ
    !Wang and Hohman
    dt_wh = 0.15d0 * sqrt(MU0*sigmin)  !myuの値max min どちらに合わせる？
    if(dt>dt4) then
    write(*,*) '****************courant 条件確認!! ********************************'
    endif
    write(*,*) '# dt mittet :', (1.0d0/sqrt(3.0d0)) * dx / cmax
    write(*,*) '# courant 2nd order scheme with taylor operator dt2', dt2!,&
    !dx/(3.0d0*0.5d0)/cmin
    write(*,*) '# courant 4th order scheme with taylor operator dt4', dt4!, &
    !(2.0d0*dx)/((3.0d0**0.5d0)*pi*cmin)
    write(*,*) '# courant dt (courant*6.0d0/7.0d0*0.999d0)', dt_ideal!,&
!     1.0d0/cmin/sqrt(1.0d0/dx**2.0d0 + 1.0d0/dy**2.0d0 + 1.0d0/dz**2.0d0)*6.0/7.0/0.999d0
!     write(*,*) '# courant Wang and Hohman dt_wh', dt_wh, 0.15d0 * sqrt(MU0*sigmax)




!CFL limit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cfl_limit = cmax*dt/dx*(3.0d0**0.5d0)
    if(cfl_limit>1) then
        write(*,*) '***************CFL limit is violated*****************'
    endif


!Fourier limit!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fourier_limit = cmax*dt/dx*pi/2.0d0*(3.0d0**0.5d0)
    if(fourier_limit>1) then
        write(*,*) '************Fourier limit is violated*******************'
    endif


!!!the number of timestep :nstepの計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     S=0.6d0  !S =[0.5:1.0]

    Nt1=int( S*nx/sqrt( 1.0d0/dx**2 + 1.0d0/dy**2 + 1.0d0/dz**2 ) )
!     Nt1= S*dble(nx)/sqrt( 1.0d0/dble(dx)**2.0d0 + 1.0d0/dble(dy)**2.0d0 + 1.0d0/dble(dz)**2.0d0 )
    Nt2=int( S*nx/sqrt( 1.0d0/dx**2 + 1.0d0/dy**2 + 1.0d0/dz**2 ) )* sqrt(sigmax/sigmin)
    Nt3=int( S*nx*(3.0d0**0.5)*sqrt(sigmax/sigmin) )
    Nt4=int( S*nx*(3.0d0**0.5d0)*cmax/cmin )

    write(*,*) '# number of timestep Nt1,Nt2,Nt3,N4',Nt1,Nt2,Nt3,Nt4

     if(Nt1 >nstep) then
    write(*,*) "******* time step nstep is violated *******"
    endif




!最大の周波数:fmaxの確認!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! grid dispersionと同じ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     cmin = sqrt(2.0d0*omega0/MU0/1.0d0)
    fmax_w = cmin /Glim /max(dx,dy,dz)

!     write(*,*) "# fmax_w",fmax_w, fmax

     if(fmax_w <fmax) then
    write(*,*), "******* fmax is violated ****fmax > cmin /Glim /max(dx,dy,dz)"
    endif
            end subroutine confirm_parameter




!//////////////////////////////////////////////////////////////////////////
!初期 ehfield set 0
!////////////////////////////////////////////////////////////////////////
subroutine set_zero_eh(EX,EY,EZ,HX,HY,HZ)
        use const_para
        implicit none

        complex(kind(0d0)), intent(out) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
        complex(kind(0d0)), intent(out) :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Ezx1(nx,ny,nz),psi_Eyx1(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Exy1(nx,ny,nz),psi_Ezy1(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Eyz1(nx,ny,nz),psi_Exz1(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Hzx1(nx,ny,nz),psi_Hyx1(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Hxy1(nx,ny,nz),psi_Hzy1(nx,ny,nz)
!         complex(kind(0d0)), intent(out)  :: psi_Hyz1(nx,ny,nz),psi_Hxz1(nx,ny,nz)

        EX(1:nx,1:ny,1:nz) = 0.0d0
        EZ(1:nx,1:ny,1:nz) = 0.0d0
        EY(1:nx,1:ny,1:nz) = 0.0d0
        HX(1:nx,1:ny,1:nz) = 0.0d0
        HZ(1:nx,1:ny,1:nz) = 0.0d0
        HY(1:nx,1:ny,1:nz) = 0.0d0

!         psi_ezx1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_eyx1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_exy1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_ezy1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_eyz1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_exz1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hzx1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hyx1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hxy1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hzy1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hyz1(1:nx,1:ny,1:nz) = 0.0d0
!         psi_hxz1(1:nx,1:ny,1:nz) = 0.0d0
end subroutine set_zero_eh





