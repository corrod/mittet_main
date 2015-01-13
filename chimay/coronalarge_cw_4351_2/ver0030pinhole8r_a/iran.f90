program main
	implicit none
	integer ::i,sp
	real(8) :: t_res(10000)
	real(8) :: signal_res(10000)
	character(16) :: file_sp

do sp = 19,103

    write(file_sp,*) sp !
	open(1000,file='./out_residual_rev/'//trim(adjustl(file_sp))//'residual_wave_diff_rev.d')
! 	open(1000,file='./out_residual_rev/19residual_wave_diff_rev.d')
	do i=1,1000
! 	read(1000,*) t_res, signal_res
	write(*,*) t_res(i),signal_res(i)
	enddo
	close(1000)
enddo

end program main