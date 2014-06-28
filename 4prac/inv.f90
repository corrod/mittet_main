
program inv

      use f95_lapack

      implicit none
      integer,parameter :: n=2
      integer :: i,j
      integer :: ipiv(1:n)
      real(8) :: r(1:n,1:n)

      write(*,'(A)') 'Input real matrix :'
      do i=1,n
          read(*,*) r(i,:)
      end do

      ! LU分解後に逆行列を計算する
      call LA_GETRF(r,ipiv)
      call LA_GETRI(r,ipiv)
      do i=1,n
          write(*,*) (r(i,j), j=1,n)
      end do

end program inv