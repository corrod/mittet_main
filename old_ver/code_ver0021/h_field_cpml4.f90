!//////////////////////////////////////////////////////////////////////////////
! H_field CPML ver
!//////////////////////////////////////////////////////////////////////////////
subroutine h_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz) !!

	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t
    complex(kind(0d0)), intent(in)   :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout):: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)

!hxcpml4
    do k = 2,nz-2
        do j = 2,ny-2
            do i = 2,nx-1
                Hx(i,j,k) = da_x(i,j,k) * Hx(i,j,k) &
                          - db_x(i,j,k) * (( c1*Ez(i,j+1,k) - c1*Ez(i,j,k) + c2*Ez(i,j+2,k) - c2*Ez(i,j-1,k) ) / khdy(j) &
				                        -  ( c1*Ey(i,j,k+1) - c1*Ey(i,j,k) + c2*Ey(i,j,k+2) - c2*Ey(i,j,k-1) ) / khdz(k))
            enddo
        enddo
    enddo

!hycpml4
    do k = 2,nz-2
        do j = 2,ny-1
            do i = 2,nx-2
                Hy(i,j,k) = da_y(i,j,k) * Hy(i,j,k) &
                          - db_y(i,j,k) * (( c1*Ex(i,j,k+1) - c1*Ex(i,j,k) + c2*Ex(i,j,k+2) - c2*Ex(i,j,k-1) ) / khdz(k) &
				                        -  ( c1*Ez(i+1,j,k) - c1*Ez(i,j,k) + c2*Ez(i+2,j,k) - c2*Ez(i-1,j,k) ) / khdx(i))
            enddo
        enddo
    enddo

!hzcpml4
    do k = 2,nz-1
        do j = 2,ny-2
            do i = 2,nx-2
                Hz(i,j,k) = da_z(i,j,k) * Hz(i,j,k)&
                          - db_z(i,j,k) * (( c1*Ey(i+1,j,k) - c1*Ey(i,j,k) + c2*Ey(i+2,j,k) - c2*Ey(i-1,j,k) ) / khdx(i) &
								   		-  ( c1*Ex(i,j+1,k) - c1*Ex(i,j,k) + c2*Ex(i,j+2,k) - c2*Ex(i,j-1,k) ) / khdy(j))
             enddo
        enddo
    enddo
end subroutine h_field_cpml4



subroutine h_field_cpml4bp(istep,t,Ex,Ey,EZ,Hx,Hy,Hz) !!

	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t
    complex(kind(0d0)), intent(in)   :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout):: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
!+-反転
!hxcpml4
    do k = 2,nz-2
        do j = 2,ny-2
            do i = 2,nx-1
                Hx(i,j,k) = da_x(i,j,k) * Hx(i,j,k) &
                          + db_x(i,j,k) * (( c1*Ez(i,j+1,k) - c1*Ez(i,j,k) + c2*Ez(i,j+2,k) - c2*Ez(i,j-1,k) ) / khdy(j) &
				                        -  ( c1*Ey(i,j,k+1) - c1*Ey(i,j,k) + c2*Ey(i,j,k+2) - c2*Ey(i,j,k-1) ) / khdz(k))
            enddo
        enddo
    enddo

!hycpml4
    do k = 2,nz-2
        do j = 2,ny-1
            do i = 2,nx-2
                Hy(i,j,k) = da_y(i,j,k) * Hy(i,j,k) &
                          + db_y(i,j,k) * (( c1*Ex(i,j,k+1) - c1*Ex(i,j,k) + c2*Ex(i,j,k+2) - c2*Ex(i,j,k-1) ) / khdz(k) &
				                        -  ( c1*Ez(i+1,j,k) - c1*Ez(i,j,k) + c2*Ez(i+2,j,k) - c2*Ez(i-1,j,k) ) / khdx(i))
            enddo
        enddo
    enddo

!hzcpml4
    do k = 2,nz-1
        do j = 2,ny-2
            do i = 2,nx-2
                Hz(i,j,k) = da_z(i,j,k) * Hz(i,j,k)&
                          + db_z(i,j,k) * (( c1*Ey(i+1,j,k) - c1*Ey(i,j,k) + c2*Ey(i+2,j,k) - c2*Ey(i-1,j,k) ) / khdx(i) &
								   		-  ( c1*Ex(i,j+1,k) - c1*Ex(i,j,k) + c2*Ex(i,j+2,k) - c2*Ex(i,j-1,k) ) / khdy(j))
             enddo
        enddo
    enddo
end subroutine h_field_cpml4bp






!     real(8)             :: CHXLY(1:nx,1:ny,1:nz), CHYLZ(1:nx,1:ny,1:nz), CHZLX(1:nx,1:ny,1:nz)
!     real(8)             :: CHXLZ(1:nx,1:ny,1:nz), CHYLX(1:nx,1:ny,1:nz), CHZLY(1:nx,1:ny,1:nz)
! !Hx4
!     do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!               CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
!               CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
!           enddo
!       enddo
!     enddo

! !Hy4
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!              CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
!              CHYLX(i,j,k) = dt / myu(i,j,k) / dx
!           enddo
!       enddo
!     enddo
! !Hz4
!     do k = 1,nz
!        do j = 1,ny
!          do i = 1,nx
!             CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
!             CHZLY(i,j,k) = dt / myu(i,j,k) / dy
!         enddo
!       enddo
!     enddo
