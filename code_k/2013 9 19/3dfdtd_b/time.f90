!**********************************************************************
! Calculation time
! 言語: FORTRAN90
! コンパイル環境: DIGITAL Visual Fortran
!
! [使用法]
! call time_start    ! 時間カウント開始場所
! call time_stop     ! 時間カウント終了場所
!**********************************************************************

module cal_time
   character(64),parameter :: time_file_name="cal_time.txt"
   real(8) :: rtc_start,rtc_stop
end module

!----------------------------------------------------------------------
! 計算開始時間
!----------------------------------------------------------------------
subroutine time_start
   use dfport    ! rtc
   use cal_time
   implicit none

   character(len=12) real_clock(3)
   integer :: date_time_start(8),date_time_stop(8)
   integer :: fp=10

   open(fp,file=time_file_name,access='SEQUENTIAL')

   call date_and_time(real_clock (1),real_clock (2),real_clock (3), &
                      date_time_start)

   write(*,"('START TIME: ',i4.4,a1,i2.2,a1,i2.2,' ', &
             i2.2,':',i2.2,' ',i2.2,'sec ',i3.3,'msec')") &
      date_time_start(1),'/',date_time_start(2),'/',date_time_start(3), &
      date_time_start(5),date_time_start(6),date_time_start(7),date_time_start(8)
   write(fp,"('START TIME: ',i4.4,a1,i2.2,a1,i2.2,' ', &
             i2.2,':',i2.2,' ',i2.2,'sec ',i3.3,'msec')") &
      date_time_start(1),'/',date_time_start(2),'/',date_time_start(3), &
      date_time_start(5),date_time_start(6),date_time_start(7),date_time_start(8)

   ! the number of seconds elapsed since 00:00:00 Greenwich mean time, January 1, 1970.
   rtc_start=rtc()

   close(fp)

   return
end subroutine

!----------------------------------------------------------------------
! 計算終了時間＆計算時間
!----------------------------------------------------------------------
subroutine time_stop
   use dfport    ! rtc
   use cal_time
   implicit none

   character(len=12) real_clock(3)
   integer :: date_time_start(8),date_time_stop(8)
   integer :: fp=10

   open(fp,file=time_file_name,access='APPEND')

   call date_and_time(real_clock (1),real_clock (2),real_clock (3), &
                      date_time_stop)

   write(*,"('STOP TIME: ',i4.4,a1,i2.2,a1,i2.2,' ', &
             i2.2,':',i2.2,' ',i2.2,'sec ',i3.3,'msec')") &
      date_time_stop(1),'/',date_time_stop(2),'/',date_time_stop(3), &
      date_time_stop(5),date_time_stop(6),date_time_stop(7),date_time_stop(8)
   write(fp,"('STOP TIME: ',i4.4,a1,i2.2,a1,i2.2,' ', &
             i2.2,':',i2.2,' ',i2.2,'sec ',i3.3,'msec')") &
      date_time_stop(1),'/',date_time_stop(2),'/',date_time_stop(3), &
      date_time_stop(5),date_time_stop(6),date_time_stop(7),date_time_stop(8)

   ! the number of seconds elapsed since 00:00:00 Greenwich mean time, January 1, 1970.
   rtc_stop=rtc()

   close(fp)

   call elapsed_time

   return
end subroutine

!----------------------------------------------------------------------
! 計算時間
! 計算時間が一日以内
!----------------------------------------------------------------------
subroutine elapsed_time
   use cal_time
   implicit none

   real(8) :: time_spent
   integer :: h,m,s
   integer :: fp=10

   open(fp,file=time_file_name,access='APPEND')

   time_spent=rtc_stop-rtc_start

   s=mod(int(time_spent),60)
   m=mod(int(time_spent/60),60)
   h=time_spent/60/60
   write(*,"('ELAPSED TIME (HOUR:MINUTE:SECOND): ',i2.2,':',i2.2,':',i2.2)") h,m,s
   write(fp,"('ELAPSED TIME (HOUR:MINUTE:SECOND): ',i2.2,':',i2.2,':',i2.2)") h,m,s

   close(fp)

   return
end subroutine
