
!****************************************************************************
! メインルーチン
!****************************************************************************
program main
   use fdtd
   implicit none

   integer :: i,j

   open(1,file='test.txt')

   call init                                  ! 初期化
   call modeling                              ! 物体のモデリング
   !call test_modeling                         ! モデリングの確認
!   call time_start                            ! 時間測定開始

   time=0.0d0

   do i=1,ntime_start
      call renew
   end do

 !  call time_stop                             ! 時間測定終了

   do i=ntime_start,ntime_end,ntime_step
      do j=1,ntime_step
         call renew                           ! 電磁界の更新処理
      end do

      call output                             ! ファイルに出力
 
   end do

   close(1)

   stop
end program

!----------------------------------------------------------------------------
! 更新処理
!----------------------------------------------------------------------------
subroutine renew
   use fdtd
   implicit none

   !call plane_wave_source                     ! 平面波による励振

   call electric_field                        ! 電界を更新
   call ecur_infinitesimal_dipole_source      ! 電流源による励振
   !call electric_boundary_condition           ! 電界の境界条件
   call absorbing_boundary_condition_for_e    ! 電界に対する吸収境界条件
   time=time+dt/2.0d0

   call magnetic_field                        ! 磁界を更新
   !call magnetic_boundary_condition           ! 磁界の境界条件
   !call absorbing_boundary_condition_for_h    ! 磁界に対する吸収境界条件
   time=time+dt/2.0d0

   return
end subroutine

!
! End of file
!
