!//////////////////////////////////////////////////////////////////////////////
! E_field CPML ver
!//////////////////////////////////////////////////////////////////////////////
subroutine media_coeff
	use const_para
	ca_x(i,j,k) = (1.0d0 - ((esigx(i)*dt)/(2.0d0*epsi)))
	implicit none
end subroutine media_coeff


subroutine e_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz,sig)
	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: sig(1:nx,1:ny,1:nz)
!     real(8)             :: etaxx(nx,ny,nz),etayy(nx,ny,nz),etazz(nx,ny,nz)
!     real(8)             :: CEXLY(1:nx,1:ny,1:nz),CEYLZ(1:nx,1:ny,1:nz),CEZLX(1:nx,1:ny,1:nz)
!     real(8)             :: CEXLZ(1:nx,1:ny,1:nz),CEYLX(1:nx,1:ny,1:nz),CEZLY(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in)    :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

    !ex_cpml
     do k = 3,nz-1
        do j = 3,ny-1
            do i = 1,nx-1
                Ex(i,j,k) = ca_x(i,j,k) * Ex(i,j,k) &
                          + cb_x(i,j,k) * (( c1*Hz(i,j,k) - c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k) - c2*Hz(i,j-2,k) ) /kedy(j) &
                          				-  ( c1*Hy(i,j,k) - c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1) - c2*Hy(i,j,k-2) ) /kedz(k))
            enddo
        enddo
    enddo

    !ey_cpml
	do k = 3,nz-1
        do j = 1,ny-1
            do i = 3,nx-1
                Ey(i,j,k) = ca_y(i,j,k) * Ey(i,j,k) &
                          + cb_y(i,j,k) * (( c1*Hx(i,j,k) - c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1) - c2*Hx(i,j,k-2) ) / kedz(k) &
                        				-  ( c1*Hz(i,j,k) - c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k) - c2*Hz(i-2,j,k) ) / kedx(i))
            enddo
        enddo
    enddo

    !ez_cpml
    do k = 1,nz-1
        do j = 3,ny-1
            do i = 3,nx-1
                Ez(i,j,k) = ca_z(i,j,k) * Ez(i,j,k) &
                          + cb_z(i,j,k) * (( c1*Hy(i,j,k) - c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k) - c2*Hy(i-2,j,k) ) / kedx(i) &
				                        -  ( c1*Hx(i,j,k) - c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k) - c2*Hx(i,j-2,k) ) / kedy(j))
            enddo
        enddo
    enddo
end subroutine e_field_cpml4


subroutine e_field_cpml4bp(istep,t,Ex,Ey,EZ,Hx,Hy,Hz,sig)
	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), intent(in) :: sig(1:nx,1:ny,1:nz)
!     real(8)             :: etaxx(nx,ny,nz),etayy(nx,ny,nz),etazz(nx,ny,nz)
!     real(8)             :: CEXLY(1:nx,1:ny,1:nz),CEYLZ(1:nx,1:ny,1:nz),CEZLX(1:nx,1:ny,1:nz)
!     real(8)             :: CEXLZ(1:nx,1:ny,1:nz),CEYLX(1:nx,1:ny,1:nz),CEZLY(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in)    :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
! +-反転
    !ex_cpml
     do k = 3,nz-1
        do j = 3,ny-1
            do i = 1,nx-1
                Ex(i,j,k) = ca_x(i,j,k) * Ex(i,j,k) &
                          - cb_x(i,j,k) * (( c1*Hz(i,j,k) - c1*Hz(i,j-1,k) + c2*Hz(i,j+1,k) - c2*Hz(i,j-2,k) ) /kedy(j) &
                          				-  ( c1*Hy(i,j,k) - c1*Hy(i,j,k-1) + c2*Hy(i,j,k+1) - c2*Hy(i,j,k-2) ) /kedz(k))
            enddo
        enddo
    enddo

    !ey_cpml
	do k = 3,nz-1
        do j = 1,ny-1
            do i = 3,nx-1
                Ey(i,j,k) = ca_y(i,j,k) * Ey(i,j,k) &
                          - cb_y(i,j,k) * (( c1*Hx(i,j,k) - c1*Hx(i,j,k-1) + c2*Hx(i,j,k+1) - c2*Hx(i,j,k-2) ) / kedz(k) &
                        				-  ( c1*Hz(i,j,k) - c1*Hz(i-1,j,k) + c2*Hz(i+1,j,k) - c2*Hz(i-2,j,k) ) / kedx(i))
            enddo
        enddo
    enddo

    !ez_cpml
    do k = 1,nz-1
        do j = 3,ny-1
            do i = 3,nx-1
                Ez(i,j,k) = ca_z(i,j,k) * Ez(i,j,k) &
                          - cb_z(i,j,k) * (( c1*Hy(i,j,k) - c1*Hy(i-1,j,k) + c2*Hy(i+1,j,k) - c2*Hy(i-2,j,k) ) / kedx(i) &
				                        -  ( c1*Hx(i,j,k) - c1*Hx(i,j-1,k) + c2*Hx(i,j+1,k) - c2*Hx(i,j-2,k) ) / kedy(j))
            enddo
        enddo
    enddo
end subroutine e_field_cpml4bp









!     !Ex4
!     do k = 1,nz
!          do j = 1,ny
!               do i = 1,nx
!                      etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)
!                      CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
!                      CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
!               enddo
!          enddo
!     enddo

!     !Ey4
!     do k = 1,nz
!         do j = 1,ny
!            do i = 1,nx
!                   etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)
!                   CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
!                   CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
!            enddo
!         enddo
!     enddo
!     !Ez4
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                    etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)
!                    CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
!                    CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
!             enddo
!         enddo
!     enddo

