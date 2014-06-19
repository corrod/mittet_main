!------------------------------------------------------------------------------
!メインルーチン
!------------------------------------------------------------------------------

program main 
  use fdtd
  implicit none

  integer :: t1,t2,t_rate,t_max,diff

  open(1,file='test.txt')

  call init             ! 初期化
  call modeling         ! 物体のモデリング
  call test_modeling   ! モデリングの確認

  call system_clock(t1)

  time=0.0d0

  print *, "dt=",dt
     do n=1,ntime_step
        print *, n
        call renew      ! 電磁界の更新処理
        call output2
        write(10,*)
        write(20,*)
        write(30,*)
     end do

  call output1

  call system_clock(t2,t_rate,t_max)
  if (t2<t1) then
     diff=t_max-t1+t2
  else
     diff=t2-t1
  end if
  print "(A,F10.3)","time it took was:", diff/dble(t_rate)

  close(1)
  close(40)
  
  stop
end program

!-------------------------------------------------------------------------------
! 更新処理
!-------------------------------------------------------------------------------

subroutine renew
  use fdtd
  implicit none

  ! call plane_wave_source                  ! 平面波による励振

  call electric_field                       ! 電界を更新
  call ecur_infinitesimal_dipole_source     ! 電流源による励振
!  call ecur_line_source
!  call electric_boundary_condition        ! 電界の境界条件
  call absorbing_boundary_condition_for_e  ! 電界に対する吸収境界条件
  time=time+dt/2.0d0

  call magnetic_field                       ! 磁界を更新
!  call magnetic_boundary_condition        ! 磁界の境界条件
!  call absorbing_boundary_condition_for_h ! 磁界に対する吸収境界条件
  time=time+dt/2.0d0

  return
end subroutine
