!!!dt,dx,dy,dzの設定cmax,cminの計算************************************************
subroutine set_d_txyz(cmax)
    use const_para
    implicit none

    integer :: Nt
    real(8) :: S !タイムステップ数を求める際の係数
    real(8) :: cwa, cfe, cair
    real(8) :: courant 
    real(8) :: cmax
    real(8) :: cmin
    real(8) :: t_max
    real(8) :: dt_max
    real(8) :: fmax_w !最大周波数（上限）
    real(8) :: dt_ideal !クーラン条件を満たすdt


    !媒質中の伝播速度計算
    cwa  = sqrt(2.0d0*omega0/myuwa/sigmawa)   !!!!cmin=sqrt(2.0d0*omega0/myuwa/maxv)
    cfe  = sqrt(2.0d0*omega0/myufe/sigmafe)   !!!!cmax=sqrt(2.0d0*omega0/myufe/maxv)
    cair = 1.0d0 / sqrt(myuair*epsiair)
    
    write(*,'(a,2e12.4)') '伝播速度cwa,cfe', cwa,cfe !海水の伝播速度出力
    write(*,'(a,e12.4)')  '伝播距離c*t'    , cwa*dt*nstep

    cmax = cwa  !max(cwa,cfe,cair)!   最大伝播速度cmax計算
  !  write(*,*) 'cmax=',cmax


    t_max = dt*nstep !   計測時間t_max
    write(*,'(a,e12.4)')    '計測時間t_max',t_max
    dt_max = (2.0d0*dx)/((3**0.5)*3.14159d0*cwa) 
    write(*,'(a,e12.4)')    '最大タイムステップ長dt_max', dt_max 
    write(*,'(a,e12.4)')    '１辺の長さm',dx*nx
    write(*,'(a,e12.4)')    'タイムステップ長dt',dt 
    write(*,'(a,3e12.4)')   'グリッドサイズdx,dy,dz',dx,dy,dz 
    write(*,'(a,e12.4)')    'tau0',tau0 !ソースの継続時間
    write(*,'(a,3i5)')      'グリッド数nx,ny,nz',nx,ny,nz 
    write(*,'(a,i5)')       '総タイムステップnstep',nstep
    write(*,'(a,e12.4)')    'omega0',omega0
    write(*,*)              'sigmawa,myuwa',sigmawa,myuwa


!!!グリッド間隔dxの設定;グリッド分散

!!!タイムステップdtの設定;クーラン条件
  !クーラン条件
    courant = 1.0d0/cmax/sqrt(1.0d0/dx**2.0d0 +1.0d0/dy**2.0d0 +1.0d0/dz**2.0d0)
    dt_ideal = courant*6.0d0/7.0d0*0.999d0

!!!総タイムステップnstepの確
     S=0.6d0

     Nt=int(S*nx/sqrt( 1.0d0/dx**2 +1.0d0/dy**2 + 1.0d0/dz**2 ) )
    write(*,*) '必要なタイムステップ数Nt   =',Nt

     if(Nt >nstep) then
    write(*,*) "******* time step nstep is violated *******\n"
    endif


!最大の周波数の確認
    cmin = sqrt(2.0d0*omega0/MU0/1.0d0)
    fmax_w = cmin /Glim /max(dx,dy,dz)

    write(*,*) "上限の周波数fmax_w",fmax_w    
     if(fmax_w <fmax) then
    write(*,*), "******* fmax is violated *******\n"
    endif
            end subroutine set_d_txyz