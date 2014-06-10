!!!dt,dx,dy,dzの設定cmax,cminの計算************************************************
subroutine set_d_txyz
    use const_para
    implicit none

    real(8) :: cwa, cfe, cair
    real(8) :: cmax
    real(8) :: t_max
    real(8) :: dt_max

    !媒質中の伝播速度計算
    cwa = sqrt(2.0d0*omega0/myuwa/sigmawa)
    cfe = sqrt(2.0d0*omega0/myufe/sigmafe)
    cair = 1.0d0 / sqrt(myuair*epsiair)
    write(*,'(a,2e12.4)') '伝播速度cwa,cfe', cwa,cfe !海水の伝播速度出力
    write(*,'(a,e12.4)')  '伝播距離c*t'    , cwa*dt*nstep

  !  cmax = max(cwa,cfe,cair)!   最大伝播速度cmax計算
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
            end subroutine set_d_txyz
