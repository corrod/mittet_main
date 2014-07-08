module para
	implicit none
	integer :: i,j,k
	integer,parameter :: nx=100,ny=100,nz=100
	real(8) :: ex(nx,ny,nz),ey(nx,ny,nz),ez(nx,ny,nz)
	real(8) :: aaa(nx),bbb(ny),ccc(nz)
end module para


subroutine abc
	use para
	implicit none

	do k=1,nz
		do j=1,ny
			do i=1,nx
	ex(i,j,k) =1000.0d0
	ey(i,j,k) =0.0d0
	ez(i,j,k) = nx *3.0d0
	enddo
	enddo
	enddo
end subroutine abc


	subroutine def
		use para
		implicit none
		do i=1,nx
			aaa(i)=0.0d0
			bbb(i)=0.0d0
			ccc(i)=0.0d0
		aaa(10) = ex(2,2,2) * 3.0d0
		bbb(2) = ey(3,3,3) *4.0d0
		ccc(5) = ez(100,100,100)
		enddo
	end subroutine def

program subr
	use para
	implicit none

	call abc
	call def
	do i=1,100
	write(*,*) aaa(i),bbb(i),ccc(i)
	enddo
end program subr
