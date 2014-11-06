subroutine source_posi

	use const_para
	implicit none
! 	read(*,*) x_source, y_source
	!送信位置 左
	x_source  = 30!x0
	y_source  = 30!y0
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
