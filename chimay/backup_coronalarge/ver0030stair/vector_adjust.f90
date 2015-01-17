subroutine vector_ajust
	use const_para
	implicit none
	integer :: nd,ios
	real(8), allocatable :: in1_1(:),in1_2(:),in1_3(:)
	real(8), allocatable :: in2_1(:),in2_2(:),in2_3(:)
	real(8), allocatable :: in3_1(:),in3_2(:),in3_3(:)

    open(51,file='aaa',action='read')
      nd=0
	    do
	        read(51,'(f12.0)',iostat=ios)
	        if (ios<0) exit
	         nd=nd+1
	    enddo
    close(51)

    allocate(in1_1(nd), in1_2(nd), in1_3(nd))
    allocate(in2_1(nd), in2_2(nd), in2_3(nd))
    allocate(in3_1(nd), in3_2(nd), in3_3(nd))


    open(52,file=trim(adjustl(file_sp))//'pattern2_1.d')
        do i=1,nd
            read(52,*) in1_1(i),in1_2(i),in1_3(i)
        enddo
    close(52)

    open(53,file=trim(adjustl(file_sp))//'pattern2_2.d')
        do i=1,nd
            read(53,*) in2_1(i),in2_2(i),in2_3(i)
        enddo
    close(53)

    open(54,file=trim(adjustl(file_sp))//'pattern2_3.d')
        do i=1,nd
            read(54,*) in3_1(i),in3_2(i),in3_3(i)
        enddo
    close(54)



    open(55,file=trim(adjustl(file_sp))//'pattern2_12diff.d')
        do i=1,nd
            write(55,*) in1_1(i), in1_2(i) - in2_2(i)
        enddo
    close(55)

    open(56,file=trim(adjustl(file_sp))//'pattern2_23diff.d')
        do i=1,nd
            write(56,*) in2_1(i), in2_2(i) - in3_2(i)
        enddo
    close(56)

    deallocate(in1_1,in1_2,in1_3)
    deallocate(in2_1,in2_2,in2_3)
    deallocate(in3_1,in3_2,in3_3)
end subroutine vector_ajust
