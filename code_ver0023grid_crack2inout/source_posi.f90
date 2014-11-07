!/////////////////////////////////////////////////////
! 送信源位置設定 source_posi.f90
! x_source を決めるとx_source2、受信点1、2、3も決まる
! L_hori x_source, x_source2間距離(横)
! L_sr x_source, recever間距離(縦)
!/////////////////////////////////////////////////////

subroutine source_posi
	use const_para
	implicit none

	!送信位置 左
! 	write(*,*) "x_source位置 :"
! 	read(*,*) x_source

 	x_source  = sp
	y_source = y0
    z_source  = plate + offset
	!送信位置 右
    x_source2 = x_source + L_hori
    y_source2 = y_source
    z_source2 = z_source
    !受信位置
    !①
    xx1 = x_source
    yy1 = y_source
    zz1 = z_source + L_sr
    !②
    xx2 = xx1
    yy2 = yy1
    zz2 = z_source - L_sr
    !③
    xx3 = xx1 + L_hori
    yy3 = yy1
    zz3 = z_source - L_sr
    !①
    !● ●
    !② ③
	end subroutine source_posi
