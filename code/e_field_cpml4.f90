!////////////////////////////////////////////////////////////////////
! CPMl用の伝播計算の係数
!///////////////////////////////////////////////////////////////////
subroutine media_coeff
	use const_para
	implicit none
		sigmax2 = 1.0d0!sig(1)
		gradmax = 1.0d0!grad(1)
		scaler = 0.01d0 * sigmax2/gradmax
		do k=1,nz
			do j=1,ny
				do i=1,nx
					sig2(i,j,k) = sig(i,j,k) - scaler * grad(i,j,k)
				enddo
			enddo
		enddo
!  //scaler = 0.01f * sigx2;

!CPML coefficient
	do k=1,nz
		do j=1,ny
			do i=1,nx
				eps2(i,j,k) = sig2(i,j,k) / 2.0d0 / omega0
				!CPML coefficient
				ca_x(i,j,k) = (1.0d0 - ((esig_x(i)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((esig_x(i)*dt)/(2.0d0*eps2(i,j,k))))
				ca_y(i,j,k) = (1.0d0 - ((esig_y(j)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((esig_y(j)*dt)/(2.0d0*eps2(i,j,k))))
				ca_z(i,j,k) = (1.0d0 - ((esig_z(k)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((esig_z(k)*dt)/(2.0d0*eps2(i,j,k))))

				da_x(i,j,k) = (1.0d0 - ((msig_x(i)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((msig_x(i)*dt)/(2.0d0*eps2(i,j,k))))
				da_y(i,j,k) = (1.0d0 - ((msig_y(j)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((msig_y(j)*dt)/(2.0d0*eps2(i,j,k))))
				da_z(i,j,k) = (1.0d0 - ((msig_z(k)*dt)/(2.0d0*eps2(i,j,k)))) &
							/ (1.0d0 + ((msig_z(k)*dt)/(2.0d0*eps2(i,j,k))))

				cb_x(i,j,k) = dt/eps2(i,j,k) /(1.0d0+(esig_x(i)*dt)/(2.0d0*eps2(i,j,k)))
				cb_y(i,j,k) = dt/eps2(i,j,k) /(1.0d0+(esig_y(j)*dt)/(2.0d0*eps2(i,j,k)))
				cb_z(i,j,k) = dt/eps2(i,j,k) /(1.0d0+(esig_z(k)*dt)/(2.0d0*eps2(i,j,k)))

				db_x(i,j,k) = dt/MU0 /(1.0d0+(msig_x(i)*dt)/(2.0d0*eps2(i,j,k)))
				db_y(i,j,k) = dt/MU0 /(1.0d0+(msig_y(j)*dt)/(2.0d0*eps2(i,j,k)))
				db_z(i,j,k) = dt/MU0 /(1.0d0+(msig_z(k)*dt)/(2.0d0*eps2(i,j,k)))

	    enddo
	  enddo
	enddo

do k=1,nz
  do j=1,ny
    do i=1,nx
    !imamu system
    ca_x(i,j,k) = (1.0d0-esig_x(i)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+esig_x(i)*dt/(2.0d0*epsi(i,j,k)))
    ca_y(i,j,k) = (1.0d0-esig_y(j)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+esig_y(j)*dt/(2.0d0*epsi(i,j,k)))
    ca_z(i,j,k) = (1.0d0-esig_z(k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+esig_z(k)*dt/(2.0d0*epsi(i,j,k)))
    cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+esig_x(i)*dt/(2.0d0*epsi(i,j,k)))
    cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+esig_y(j)*dt/(2.0d0*epsi(i,j,k)))
    cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+esig_z(k)*dt/(2.0d0*epsi(i,j,k)))

    !saito system
   ! ca_x(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! ca_y(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! ca_z(i,j,k) = (1.0d0-sig(i,j,k)*dt/(2.0d0*epsi(i,j,k))) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_x(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_y(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   ! cb_z(i,j,k) = (dt/epsi(i,j,k)) / (1.0d0+sig(i,j,k)*dt/(2.0d0*epsi(i,j,k)))
   enddo
  enddo
enddo

   do k=1,nz
    do j=1,ny
        do i=1,nx
        !imamu system
        da_x(i,j,k) = (1.0d0 - ((msig_x(i)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msig_x(i)*dt)/(2.0d0*epsi(i,j,k))))
        da_y(i,j,k) = (1.0d0 - ((msig_y(j)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msig_y(j)*dt)/(2.0d0*epsi(i,j,k))))
        da_z(i,j,k) = (1.0d0 - ((msig_z(k)*dt)/(2.0d0*epsi(i,j,k)))) / (1.0d0 + ((msig_z(k)*dt)/(2.0d0*epsi(i,j,k))))
        !!!****imamu systemではmyu→MU0 myu→eps2
        db_x(i,j,k) = dt/MU0/(1.0d0+(msig_x(i)*dt)/(2.0d0*epsi(i,j,k)))
        db_y(i,j,k) = dt/MU0/(1.0d0+(msig_y(j)*dt)/(2.0d0*epsi(i,j,k)))
        db_z(i,j,k) = dt/MU0/(1.0d0+(msig_z(k)*dt)/(2.0d0*epsi(i,j,k)))

         !saito system
    !   da_x(i,j,k) = (1.0d0-(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k))) !sig=σ*
    !   da_y(i,j,k) = (1.0d0-(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k)))!導磁率σ
    !   da_z(i,j,k) = (1.0d0-(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k))) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   saito system
    !   db_x(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   db_y(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
    !   db_z(i,j,k) = (dt/myu(i,j,k)) / (1.0d0+(sig(i,j,k)*dt)/(2.0d0*myu(i,j,k)))
!!!sig_x(i)?sig(i,j,k)?
	        enddo
	    enddo
	enddo

end subroutine media_coeff




!//////////////////////////////////////////////////////////////////////////////
! E_field CPML ver
!//////////////////////////////////////////////////////////////////////////////
subroutine e_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz)
! subroutine e_field_cpml4(istep,t,Ex,Ey,EZ,Hx,Hy,Hz,sig)
	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!     real(8), intent(in) :: sig(1:nx,1:ny,1:nz)
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



subroutine e_field_cpml4bp(istep,t,Ex,Ey,EZ,Hx,Hy,Hz) !!
! subroutine e_field_cpml4bp(istep,t,Ex,Ey,EZ,Hx,Hy,Hz,sig) !!
	use const_para
    implicit none

    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!     real(8), intent(in) :: sig(1:nx,1:ny,1:nz)  !!
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





!     real(8)             :: etaxx(nx,ny,nz),etayy(nx,ny,nz),etazz(nx,ny,nz)
!     real(8)             :: CEXLY(1:nx,1:ny,1:nz),CEYLZ(1:nx,1:ny,1:nz),CEZLX(1:nx,1:ny,1:nz)
!     real(8)             :: CEXLZ(1:nx,1:ny,1:nz),CEYLX(1:nx,1:ny,1:nz),CEZLY(1:nx,1:ny,1:nz)





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

