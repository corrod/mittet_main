subroutine cerjan_e(ex,ey,ez)
	use const_para
	implicit none
		integer, parameter :: beta = 20
		real(8), parameter :: alpha = 5.5d-2
	    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)

	    do k=1,nz
	    	do j=1,ny
	    		do i=1,beta
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=1,beta
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,beta
	    	do j=1,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    		enddo	
	    	enddo
	    enddo



	    do k=1,nz
	    	do j=1,ny
	    		do i=nx-beta,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=ny-beta,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=nz-beta,nz
	    	do j=1,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

end subroutine cerjan_e










subroutine cerjan_h(hx,hy,hz)
	use const_para
	implicit none
		integer, parameter :: beta = 20
		real(8), parameter :: alpha = 5.5d-2
	    complex(kind(0d0)), intent(inout) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

	    do k=1,nz
	    	do j=1,ny
	    		do i=1,beta
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=1,beta
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,beta
	    	do j=1,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    		enddo	
	    	enddo
	    enddo



	    do k=1,nz
	    	do j=1,ny
	    		do i=nx-beta,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=ny-beta,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=nz-beta,nz
	    	do j=1,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

end subroutine cerjan_h











subroutine cerjan2_e(ex,ey,ez)
	use const_para
	implicit none
		integer, parameter :: beta = 20
		real(8), parameter :: alpha = 1.0d-4  !ok!1.0d-4 maybe
	    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)

	    do k=1,nz
	    	do j=1,ny
	    		do i=1,beta
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=1,beta
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,beta
	    	do j=1,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    		enddo	
	    	enddo
	    enddo



	    do k=1,nz
	    	do j=1,ny
	    		do i=nx-beta,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=ny-beta,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=nz-beta,nz
	    	do j=1,ny
	    		do i=1,nx
	    			ex(i,j,k)= ex(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			ey(i,j,k)= ey(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			ez(i,j,k)= ez(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

end subroutine cerjan2_e



















subroutine cerjan2_h(hx,hy,hz)
	use const_para
	implicit none
		integer, parameter :: beta = 20
		real(8), parameter :: alpha =1.0d-4 !ok!1.0d-3 maybe
	    complex(kind(0d0)), intent(inout) :: hx(nx,ny,nz),hy(nx,ny,nz),hz(nx,ny,nz)

	    do k=1,nz
	    	do j=1,ny
	    		do i=1,beta
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-i)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=1,beta
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-j)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,beta
	    	do j=1,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-k)**2.0d0)
	    		enddo	
	    	enddo
	    enddo



	    do k=1,nz
	    	do j=1,ny
	    		do i=nx-beta,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(nx-i)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=1,nz
	    	do j=ny-beta,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(ny-j)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

	    do k=nz-beta,nz
	    	do j=1,ny
	    		do i=1,nx
	    			hx(i,j,k)= hx(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			hy(i,j,k)= hy(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    			hz(i,j,k)= hz(i,j,k)*exp(-alpha*dble(beta-(nz-k)+1)**2.0d0)
	    		enddo	
	    	enddo
	    enddo

end subroutine cerjan2_h