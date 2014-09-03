program main
    implicit none
    character(16) :: a, b
    open(11,file='../code_ver0021/out1/prac_out.dat')
    read(11,*) a,b
    close(11)

    write(*,*) a ,b
    write(*,*) b

end program main
