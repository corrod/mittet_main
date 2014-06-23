!!!仮想領域での電磁場の計算E-field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  subroutine の統合
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Efield(istep,t,Je,Ex,Ey,EZ,Hx,Hy,Hz,sig)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
    real(8), parameter  :: c1 = 1.125d0, c2 = -0.04167d0 ! from Tayor expansion
    real(8), intent(in) :: Je(nstep)
    real(8), intent(in) :: sig(1:nx,1:ny,1:nz)
    real(8)             :: etaxx(nx,ny,nz),etayy(nx,ny,nz),etazz(nx,ny,nz)
    real(8)             :: CEXLY(1:nx,1:ny,1:nz),CEYLZ(1:nx,1:ny,1:nz),CEZLX(1:nx,1:ny,1:nz)
    real(8)             :: CEXLZ(1:nx,1:ny,1:nz),CEYLX(1:nx,1:ny,1:nz),CEZLY(1:nx,1:ny,1:nz)
!     complex(kind(0d0)), intent(inout):: Ex(-1:nx+2,-1:ny+2,-1:nz+2),Ey(-1:nx+2,-1:ny+2,-1:nz+2),Ez(-1:nx+2,-1:ny+2,-1:nz+2)
!     complex(kind(0d0)), intent(in)   :: Hx(-1:nx+2,-1:ny+2,-1:nz+2),Hy(-1:nx+2,-1:ny+2,-1:nz+2),Hz(-1:nx+2,-1:ny+2,-1:nz+2)
    complex(kind(0d0)), intent(inout):: Ex(nx,ny,nz),Ey(nx,ny,nz),Ez(nx,ny,nz)
    complex(kind(0d0)), intent(in)   :: Hx(nx,ny,nz),Hy(nx,ny,nz),Hz(nx,ny,nz)
                                      !   complex(kind(0d0)), intent(in) :: sigxx(nx,ny,nz)
                                      !     real(8)             :: alpha(ln,ln)
                                      !     real(8)             :: temp1(0:ln),temp2(0:ln)
                                       !   temp1(0) = 0.0d0
                                        !  temp2(0) = 0.0d0

                                          !Holberg optimization scheme
                                       !   alpha(1,1)    = 1.00235d0
                                          !alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
                                          !alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
                                          !alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)

                                          !Taylor expansion
                                      !     alpha(1,1)   = 1.0d0
                                      !     alpha(2,1:2) = (/1.12500d0,-0.04167d0/)
                                      !     alpha(3,1:3) = (/1.17188d0,-0.06510d0,0.00469d0/)
                                      !     alpha(4,1:4) = (/1.19629d0,-0.07975d0,-0.00070d0/)


!ver1   444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
!ln=2 c1 = 1.125d0, c2 = -0.04167d0 ! from Tayor expansion
    !Ex4
    do k = 1,nz
         do j = 1,ny
              do i = 1,nx
                     etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigxx(i,j,k)
                     CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
                     CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
              enddo
         enddo
    enddo
    !Ex4波動伝播計算
    do k = 3,nz-1
        do j = 3,ny-1
            do i = 1,nx-1
                Ex(i,j,k) = Ex(i,j,k)&
                         + CEXLY(i,j,k)*( c1*Hz(i,j,k)-c1*Hz(i,j-1,k) +c2*Hz(i,j+1,k) - c2*Hz(i,j-2,k) )&
                         + CEXLZ(i,j,k)*( c1*Hy(i,j,k)-c1*Hy(i,j,k-1) +c2*Hy(i,j,k+1) - c2*Hy(i,j,k-2) )
            enddo
        enddo
    enddo

    !ソース項
    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)

    !Ey4
    do k = 1,nz
        do j = 1,ny
           do i = 1,nx
                  etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigyy(i,j,k)
                  CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
                  CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
           enddo
        enddo
    enddo
    !Ey4波動伝播計算
    do k = 3,nz-1
        do j = 1,ny-1
            do i = 3,nx-1
                Ey(i,j,k) = Ey(i,j,k)&
                         + CEYLZ(i,j,k)*( c1*Hx(i,j,k)-c1*Hx(i,j,k-1) +c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2) )&
                         + CEYLX(i,j,k)*( c1*Hz(i,j,k)-c1*Hz(i-1,j,k) +c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k) )
            enddo
        enddo
    enddo

    !ソース項
    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)


    !Ez4
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                   etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigzz(i,j,k)
                   CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
                   CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
            enddo
        enddo
    enddo
    !Ez4波動伝播計算
    do k = 1,nz-1
        do j = 3,ny-1
            do i = 3,nx-1
                Ez(i,j,k) = Ez(i,j,k)&
                          + CEZLX(i,j,k)*( c1*Hy(i,j,k) - c1*Hy(i-1,j,k) +c2*Hy(i+1,j,k) - c2*Hy(i-2,j,k) )&
                          + CEZLY(i,j,k)*( c1*Hx(i,j,k) - c1*Hx(i,j-1,k) +c2*Hx(i,j+1,k) - c2*Hx(i,j-2,k) )
            enddo
        enddo
    enddo

    !ソース項
!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)










!NORMAL version2 ln=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !    !Ex係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!         do j = 1,ny
!              do i = 1,nx
!                     etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigxx(i,j,k)
!                     CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
!                     CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
!              enddo
!         enddo
!    enddo

!    !Ex波動伝播計算
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Ex(i,j,k) = Ex(i,j,k)&
!                           + CEXLY(i,j,k)*(Hz(i,j,k)-Hz(i,j-1,k))&
!                           + CEXLZ(i,j,k)*(Hy(i,j,k)-Hy(i,j,k-1))
!            enddo
!        enddo
!    enddo

!    !ソース項
!    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)

! !    !Ey係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!           do i = 1,nx
!                  etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigyy(i,j,k)
!                  CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
!                  CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
!           enddo
!        enddo
!    enddo

!    !Ey波動伝播計算
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Ey(i,j,k) = Ey(i,j,k)&
!                           + CEYLZ(i,j,k)*(Hx(i,j,k)-Hx(i,j,k-1))&
!                           + CEYLX(i,j,k)*(Hz(i,j,k)-Hz(i-1,j,k))
!            enddo
!        enddo
!    enddo

!    !ソース項
!    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)


! !    !Ez係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                   etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigzz(i,j,k)
!                   CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
!                   CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
!            enddo
!        enddo
!    enddo

!    !Ez波動伝播計算
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                Ez(i,j,k) = Ez(i,j,k)&
!                           + CEZLX(i,j,k)*(Hy(i,j,k)-Hy(i-1,j,k))&
!                           + CEZLY(i,j,k)*(Hx(i,j,k)-Hx(i,j-1,k))
!             enddo
!        enddo
 !  enddo

!    !ソース項
!!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)


!ver2 imamu c style444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
!     !Ex4
!     do k = 1,nz
!          do j = 1,ny
!               do i = 1,nx
!                      etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigxx(i,j,k)
!                      CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
!                      CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
!               enddo
!          enddo
!     enddo
!     !Ex4波動伝播計算
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Ex(i,j,k) = Ex(i,j,k)&
!                          + CEXLY(i,j,k)*( c1*Hz(i,j,k)-c1*Hz(i,j-1,k) +c2*Hz(i,j+1,k) - c2*Hz(i,j-2,k) )&
!                          + CEXLZ(i,j,k)*( c1*Hy(i,j,k)-c1*Hy(i,j,k-1) +c2*Hy(i,j,k+1) - c2*Hy(i,j,k-2) )
!             enddo
!         enddo
!     enddo

!     !ソース項
!     !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)

!     !Ey4
!     do k = 1,nz
!         do j = 1,ny
!            do i = 1,nx
!                   etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigyy(i,j,k)
!                   CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
!                   CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
!            enddo
!         enddo
!     enddo
!     !Ey4波動伝播計算
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Ey(i,j,k) = Ey(i,j,k)&
!                          + CEYLZ(i,j,k)*( c1*Hx(i,j,k)-c1*Hx(i,j,k-1) +c2*Hx(i,j,k+1)-c2*Hx(i,j,k-2) )&
!                          + CEYLX(i,j,k)*( c1*Hz(i,j,k)-c1*Hz(i-1,j,k) +c2*Hz(i+1,j,k)-c2*Hz(i-2,j,k) )
!             enddo
!         enddo
!     enddo

!     !ソース項
!     !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)


!     !Ez4
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                    etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigzz(i,j,k)
!                    CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
!                    CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
!             enddo
!         enddo
!     enddo
!     !Ez4波動伝播計算
!     do k = 1,nz
!         do j = 1,ny
!             do i = 1,nx
!                 Ez(i,j,k) = Ez(i,j,k)&
!                           + CEZLX(i,j,k)*( c1*Hy(i,j,k) - c1*Hy(i-1,j,k) +c2*Hy(i+1,j,k) - c2*Hy(i-2,j,k) )&
!                           + CEZLY(i,j,k)*( c1*Hx(i,j,k) - c1*Hx(i,j-1,k) +c2*Hx(i,j+1,k) - c2*Hx(i,j-2,k) )
!             enddo
!         enddo
!     enddo

!     !ソース項
! !   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)






!NORMAL ln=1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Ex係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!         do j = 1,ny
!              do i = 1,nx
!                     etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigxx(i,j,k)
!                     CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
!                     CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
!              enddo
!         enddo
!    enddo
!
!    !Ex波動伝播計算
!    do k = 2,nz-1
!        do j = 2,ny-1
!            do i = 1,nx-1
!                Ex(i,j,k) = Ex(i,j,k)&
                 !          + CEXLY(i,j,k)*(Hz(i,j,k)-Hz(i,j-1,k))&
                 !          + CEXLZ(i,j,k)*(Hy(i,j,k)-Hy(i,j,k-1))
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)
!
!    !Ey係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!           do i = 1,nx
!                  etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigyy(i,j,k)
!                  CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
!                  CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
!           enddo
!        enddo
!    enddo
!
!    !Ey波動伝播計算
!    do k = 2,nz-1
!        do j = 1,ny-1
!            do i = 2,nx-1
!                Ey(i,j,k) = Ey(i,j,k)&
                     !      + CEYLZ(i,j,k)*(Hx(i,j,k)-Hx(i,j,k-1))&
                     !      + CEYLX(i,j,k)*(Hz(i,j,k)-Hz(i-1,j,k))
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)
!
!
!    !Ez係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                   etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigzz(i,j,k)
!                   CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
!                   CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
!            enddo
!        enddo
!    enddo
!
!    !Ez波動伝播計算
!    do k = 1,nz-1
!        do j = 2,ny-1
!            do i = 2,nx-1
!                Ez(i,j,k) = Ez(i,j,k)&
                    !       + CEZLX(i,j,k)*(Hy(i,j,k)-Hy(i-1,j,k))&
                    !       + CEZLY(i,j,k)*(Hx(i,j,k)-Hx(i,j-1,k))
 !            enddo
!        enddo
!    enddo
!
!    !ソース項
!!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)



!ボツ ln=2 original ver
!!!!!!! ln=1~4version !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Ex係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!         do j = 1,ny
!              do i = 1,nx
!                     etaxx(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigxx(i,j,k)
!                     CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
!                     CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
!              enddo
!         enddo
!    enddo
!    !Ex波動伝播計算
!    do k = 2,nz-1
!        do j = 2,ny-1
!            do i = 1,nx-1
!                    do l=1,ln
!                        temp1(l) = temp1(l-1) + alpha(ln,l) * (Hz(i,j+(l-1),k) - Hz(i,j-l,k))
!                        temp2(l) = temp2(l-1) + alpha(ln,l) * (Hy(i,j,k+(l-1)) - Hy(i,j,k-l))
!                    enddo
!                Ex(i,j,k) = Ex(i,j,k) + CEXLY(i,j,k)*temp1(ln) + CEXLZ(i,j,k)*temp2(ln)
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)
!
!    !Ey係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!           do i = 1,nx
!                  etayy(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigyy(i,j,k)
!                  CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
!                  CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
!           enddo
!        enddo
!    enddo
!    !Ey波動伝播計算
!    do k = 2,nz-1
!        do j = 1,ny-1
!            do i = 2,nx-1
!                    do l=1,ln
!                        temp1(l) = temp1(l-1) + alpha(ln,l) * (Hx(i,j,k+(l-1)) - Hx(i,j,k-l))
!                        temp2(l) = temp2(l-1) + alpha(ln,l) * (Hz(i+(l-1),j,k) - Hz(i-l,j,k))
!                    enddo
!                Ey(i,j,k) = Ey(i,j,k) + CEYLZ(i,j,k)*temp1(ln) + CEYLX(i,j,k)*temp2(ln)
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)
!
!
!    !Ez係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do k = 1,nz
!        do j = 1,ny
!            do i = 1,nx
!                   etazz(i,j,k) = 2.0d0 * omega0 * sig(i,j,k)!sigzz(i,j,k)
!                   CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
!                   CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
!            enddo
!        enddo
!    enddo
!    !Ez波動伝播計算
!    do k = 1,nz-1
!        do j = 2,ny-1
!            do i = 2,nx-1
!                    do l=1,ln
!                        temp1(l) = temp1(l-1) + alpha(ln,l) * (Hy(i+(l-1),j,k) - Hy(i-l,j,k))
!                        temp2(l) = temp2(l-1) + alpha(ln,l) * (Hx(i,j+(l-1),k) - Hx(i,j-l,k))
!                    enddo
!                Ez(i,j,k) = Ez(i,j,k) + CEZLX(i,j,k)*temp1(ln) + CEZLY(i,j,k)*temp2(ln)
!            enddo
!        enddo
!    enddo
!
!    !ソース項
!!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)

end subroutine Efield
