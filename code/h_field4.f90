!仮想領域での磁場計算 Hfield !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine の統合
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  Hfield(istep,t,Jh,Ex,Ey,Ez,Hx,Hy,Hz,myu)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t
    real(8), intent(in) :: myu(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(in) :: Jh(nstep)
    real(8)             :: CHXLY(1:nx,1:ny,1:nz), CHYLZ(1:nx,1:ny,1:nz), CHZLX(1:nx,1:ny,1:nz)
    real(8)             :: CHXLZ(1:nx,1:ny,1:nz), CHYLX(1:nx,1:ny,1:nz), CHZLY(1:nx,1:ny,1:nz)
    complex(kind(0d0)), intent(in)   :: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout):: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
!     complex(kind(0d0)), intent(in)   :: Ex(-1:nx+2,-1:ny+2,-1:nz+2),Ey(-1:nx+2,-1:ny+2,-1:nz+2),Ez(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)), intent(inout):: Hx(-1:nx+2,-1:ny+2,-1:nz+2),Hy(-1:nx+2,-1:ny+2,-1:nz+2),Hz(-1:nx+2,-1:ny+2,-1:nz+2)
!     real(8)             :: alpha(ln,ln)
!     real(8)             :: temp1(0:ln), temp2(0:ln)
!     temp1(0) = 0.0d0
!     temp2(0) = 0.0d0

    !Holberg optimization scheme
  !   alpha(1,1) = 1.00235d0
  !  alpha(2,1:2) = (/1.14443d0, -0.04886d0/)
  !  alpha(3,1:3) = (/1.20282d0, -0.08276d0, 0.00950d0/)
  !  alpha(4,1:4) = (/1.23041d0, -0.10313d0, 0.02005d0, -0.00331d0/)

    !Taylor expansion
!     alpha(1,1)   = 1.0d0
!     alpha(2,1:2) = (/1.12500d0,-0.04167d0/)
!     alpha(3,1:3) = (/1.17188d0,-0.06510d0,0.00469d0/)
!     alpha(4,1:4) = (/1.19629d0,-0.07975d0,-0.00070d0/)



!ver1 44444444444444444444444444444444444444444444444444444444444444444444444444444444444444
!ln=2 c1 = 1.125d0, c2 = -0.04167d0 ! from Tayor expansion
!Hx4
    do k = 1,nz
       do j = 1,ny
           do i = 1,nx
              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
          enddo
      enddo
    enddo
!Hx4波動伝播
    do k = 2,nz-2
        do j = 2,ny-2
            do i = 2,nx-1
                Hx(i,j,k) = Hx(i,j,k)&
                         + CHXLY(i,j,k)*( c1*Ez(i,j+1,k) - c1*Ez(i,j,k) + c2*Ez(i,j+2,k) - c2*Ez(i,j-1,k) )&
                         + CHXLZ(i,j,k)*( c1*Ey(i,j,k+1) - c1*Ey(i,j,k) + c2*Ey(i,j,k+2) - c2*Ey(i,j,k-1) )
            enddo
        enddo
    enddo

    !ソース項
!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)

!Hy4
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
          enddo
      enddo
    enddo
!Hy4波動伝播計算
    do k = 2,nz-2
        do j = 2,ny-1
            do i = 2,nx-2
                Hy(i,j,k) = Hy(i,j,k)&
                         + CHYLZ(i,j,k)*( c1*Ex(i,j,k+1) - c1*Ex(i,j,k) + c2*Ex(i,j,k+2) - c2*Ex(i,j,k-1) )&
                         + CHYLX(i,j,k)*( c1*Ez(i+1,j,k) - c1*Ez(i,j,k) + c2*Ez(i+2,j,k) - c2*Ez(i-1,j,k) )
            enddo
        enddo
    enddo

    !ソース項
  !   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)

!Hz4
    do k = 1,nz
       do j = 1,ny
         do i = 1,nx
            CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
            CHZLY(i,j,k) = dt / myu(i,j,k) / dy
        enddo
      enddo
    enddo
!Hz4波動伝播計算
    do k = 2,nz-1
        do j = 2,ny-2
            do i = 2,nx-2
                Hz(i,j,k) = Hz(i,j,k)&
                         + CHZLX(i,j,k)* ( c1*Ey(i+1,j,k) - c1*Ey(i,j,k) + c2*Ey(i+2,j,k) - c2*Ey(i-1,j,k) )&
                         + CHZLY(i,j,k)* ( c1*Ex(i,j+1,k) - c1*Ex(i,j,k) + c2*Ex(i,j+2,k) - c2*Ex(i,j-1,k) )
             enddo
        enddo
    enddo

    !ソース項
!     Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)





! NORMAL version2 ln=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Hx係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!           do i = 1,nx
!              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
!              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
!          enddo
!      enddo
!    enddo

!    !Hx波動伝播計
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Hx(i,j,k) = Hx(i,j,k)&
!                           + CHXLY(i,j,k)* (Ez(i,j+1,k)-Ez(i,j,k))&
!                           + CHXLZ(i,j,k)* (Ey(i,j,k+1)-Ey(i,j,k))
!            enddo
!        enddo
!    enddo

!    !ソース項
! !   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)


! !Hy係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
!             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
!          enddo
!      enddo
!    enddo

!    !Hy波動伝播計算
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Hy(i,j,k) = Hy(i,j,k)&
!                          + CHYLZ(i,j,k)*(Ex(i,j,k+1)-Ex(i,j,k))&
!                          + CHYLX(i,j,k)*(Ez(i+1,j,k)-Ez(i,j,k))
!            enddo
!        enddo
!    enddo

!    !ソース項
! !   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)


! !Hz係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!         do i = 1,nx
!            CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
!            CHZLY(i,j,k) = dt / myu(i,j,k) / dy
!        enddo
!      enddo
!    enddo

!    !Hz波動伝播計算
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Hz(i,j,k) = Hz(i,j,k)&
!                          + CHZLX(i,j,k)*(Ey(i+1,j,k)-Ey(i,j,k))&
!                          + CHZLY(i,j,k)*(Ex(i,j+1,k)-Ex(i,j,k))
!             enddo
!        enddo
!    enddo

!    !ソース項
!    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)
! !

!ver2 imamu c style44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
!     !Hx4
!     do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!               CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
!               CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
!           enddo
!       enddo
!     enddo
!     !Hx4波動伝播
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Hx(i,j,k) = Hx(i,j,k)&
!                          + CHXLY(i,j,k)*( c1*Ez(i,j+1,k) - c1*Ez(i,j,k) + c2*Ez(i,j+2,k) - c2*Ez(i,j-1,k) )&
!                          + CHXLZ(i,j,k)*( c1*Ey(i,j,k+1) - c1*Ey(i,j,k) + c2*Ey(i,j,k+2) - c2*Ey(i,j,k-1) )
!             enddo
!         enddo
!     enddo

!     !ソース項
! !   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)

!   !Hy4
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!              CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
!              CHYLX(i,j,k) = dt / myu(i,j,k) / dx
!           enddo
!       enddo
!     enddo
!     !Hy4波動伝播計算
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Hy(i,j,k) = Hy(i,j,k)&
!                          + CHYLZ(i,j,k)*( c1*Ex(i,j,k+1) - c1*Ex(i,j,k) + c2*Ex(i,j,k+2) - c2*Ex(i,j,k-1) )&
!                          + CHYLX(i,j,k)*( c1*Ez(i+1,j,k) - c1*Ez(i,j,k) + c2*Ez(i+2,j,k) - c2*Ez(i-1,j,k) )
!             enddo
!         enddo
!     enddo

!     !ソース項
!   !   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)

!   !Hz4
!     do k = 1,nz
!        do j = 1,ny
!          do i = 1,nx
!             CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
!             CHZLY(i,j,k) = dt / myu(i,j,k) / dy
!         enddo
!       enddo
!     enddo
!     !Hz4波動伝播計算
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Hz(i,j,k) = Hz(i,j,k)&
!                          + CHZLX(i,j,k)* ( c1*Ey(i+1,j,k) - c1*Ey(i,j,k) + c2*Ey(i+2,j,k) - c2*Ey(i-1,j,k) )&
!                          + CHZLY(i,j,k)* ( c1*Ex(i,j+1,k) - c1*Ex(i,j,k) + c2*Ex(i,j+2,k) - c2*Ex(i,j-1,k) )
!              enddo
!         enddo
!     enddo

!     !ソース項
!     Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)




! NORMAL ln=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Hx係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!           do i = 1,nx
!              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
!              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
!          enddo
!      enddo
!    enddo
!
!    !Hx波動伝播計
!    do k = 1,nz-1
!        do j = 1,ny-1
!            do i = 2,nx-1
!                Hx(i,j,k) = Hx(i,j,k)&
!                           + CHXLY(i,j,k)* (Ez(i,j+1,k)-Ez(i,j,k))&
!                           + CHXLZ(i,j,k)* (Ey(i,j,k+1)-Ey(i,j,k))
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)
!
!
! !Hy係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
!             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
!          enddo
!      enddo
!    enddo
!
!    !Hy波動伝播計算
!    do k = 1,nz-1
!        do j = 2,ny-1
!            do i = 1,nx-1
!                Hy(i,j,k) = Hy(i,j,k)&
                        !  + CHYLZ(i,j,k)*(Ex(i,j,k+1)-Ex(i,j,k))&
                        !  + CHYLX(i,j,k)*(Ez(i+1,j,k)-Ez(i,j,k))
!            enddo
!        enddo
!    enddo
!
!    !ソース項
! !   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)
!
!
! !Hz係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!         do i = 1,nx
!            CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
!            CHZLY(i,j,k) = dt / myu(i,j,k) / dy
!        enddo
!      enddo
!    enddo
!
!    !Hz波動伝播計算
!    do k = 2,nz-1
!        do j = 1,ny-1
!            do i = 1,nx-1
!                Hz(i,j,k) = Hz(i,j,k)&
!                          + CHZLX(i,j,k)*(Ey(i+1,j,k)-Ey(i,j,k))&
!                          + CHZLY(i,j,k)*(Ex(i,j+1,k)-Ex(i,j,k))
!             enddo
!        enddo
!    enddo
!
!    !ソース項
!    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)
!
!
!

!ボツ ln=2 オリジナルver
  !    !Hx係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!           do i = 1,nx
!              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
!              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
!          enddo
!      enddo
!    enddo
!
!    !Hx波動伝播計
!    do k = 1,nz-1
!        do j = 1,ny-1
!            do i = 2,nx-1
!                do l = 1,ln
!                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ez(i,j+l,k) - Ez(i,j-(l-1),k))
!                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ey(i,j,k+l) - Ey(i,j,k-(l-1)))
!                enddo
!                Hx(i,j,k) = Hx(i,j,k) + CHXLY(i,j,k)*temp1(ln) + CHXLZ(i,j,k)*temp2(ln)
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)
!
!
! !Hy係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
!             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
!          enddo
!      enddo
!    enddo
!
!    !Hy波動伝播計算
!    do k = 1,nz-1
!        do j = 2,ny-1
!            do i = 1,nx-1
!                do l = 1,ln
!                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ex(i,j,k+l) - Ex(i,j,k-(l-1)))
!                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ez(i+l,j,k) - Ez(i-(l-1),j,k))
!                enddo
!                Hy(i,j,k) = Hy(i,j,k) + CHYLZ(i,j,k)*temp1(ln) + CHYLX(i,j,k)*temp2(ln)
!            enddo
!        enddo
!    enddo
!
!    !ソース項
! !   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)
!
!
! !Hz係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!       do j = 1,ny
!         do i = 1,nx
!            CHZLX(i,j,k) = - dt / myu(i,j,k) / dx
!            CHZLY(i,j,k) = dt / myu(i,j,k) / dy
!        enddo
!      enddo
!    enddo

!    !Hz波動伝播計算
!    do k = 2,nz-1
!        do j = 1,ny-1
!            do i = 1,nx-1
!                do l = 1,ln
!                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ey(i+l,j,k) - Ey(i-(l-1),j,k))
!                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ex(i,j+l,k) - Ex(i,j-(l-1),k))
!                enddo
!                Hz(i,j,k) = Hz(i,j,k) + CHZLX(i,j,k)*temp1(ln) + CHZLY(i,j,k)*temp2(ln)
!             enddo
!        enddo
!    enddo

!    !ソース項
!    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)

            end subroutine Hfield
