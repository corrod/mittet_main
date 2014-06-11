!!!仮想領域での電磁場の計算*****************************************************
! operator half lengthの変更
! alpha(1),alpha(2),...の決め方は？
!電場計算*********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ex-field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EXFIELD(istep,t,Je,Ex,Hy,Hz,sigma)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)), intent(in) :: sigmaxx(nx,ny,nz)
    real(8), intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    real(8),intent(in)  :: sigma(nx,ny,nz)
    real(8) :: etaxx(nx,ny,nz)
    real(8) :: CEXLY(nx,ny,nz)
    real(8) :: CEXLZ(nx,ny,nz)
    real(8) :: alpha(ln,ln)
    real(8) :: temp1(0:ln),temp2(0:ln)
    temp1(0) = 0.0d0
    temp2(0) = 0.0d0

    !Holberg optimization scheme
    alpha(1,1) = 1.00235d0
    !alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
    !alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
    !alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)


    !係数の設定
    do k = 1,nz
         do j = 1,ny
              do i = 1,nx
                    etaxx(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmaxx(i,j,k)
                    CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy 
                    CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
              enddo
         enddo
    enddo


    !波動伝播計算
    do k = 2,nz-1
        do j = 2,ny-1
            do i = 1,nx-1
                    do l=1,ln
                        temp1(l) = temp1(l-1) + alpha(ln,l) * (Hz(i,j+(l-1),k) - Hz(i,j-l,k))
                        temp2(l) = temp2(l-1) + alpha(ln,l) * (Hy(i,j,k+(l-1)) - Hy(i,j,k-l))
                    enddo
                        Ex(i,j,k) = Ex(i,j,k) + CEXLY(i,j,k)*temp1(ln) + CEXLZ(i,j,k)*temp2(ln)
            enddo
        enddo
    enddo

    !ソース項
    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)
            end subroutine EXFIELD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ey-field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EYFIELD(istep,t,Je,Ey,Hz,Hx,sigma)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)), intent(in) :: sigmayy(nx,ny,nz)
    real(8), intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ey(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hz(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    real(8), intent(in) :: sigma(nx,ny,nz)
    real(8) :: etayy(nx,ny,nz)
    real(8) :: CEYLZ(nx,ny,nz)
    real(8) :: CEYLX(nx,ny,nz)
    real(8) :: alpha(ln,ln)
    real(8) :: temp1(0:ln),temp2(0:ln)
    temp1(0) = 0.0d0
    temp2(0) = 0.0d0

   !Holberg optimization scheme
    alpha(1,1) = 1.00235d0
    !alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
    !alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
    !alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)


    !係数の設定
    do k = 1,nz
        do j = 1,ny
           do i = 1,nx
              etayy(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmayy(i,j,k)
              CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
              CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
           enddo
        enddo
    enddo


    !波動伝播計算
    do k = 2,nz-1
        do j = 1,ny-1
            do i = 2,nx-1
                do l=1,ln
                temp1(l) = temp1(l-1) + alpha(ln,l) * (Hx(i,j,k+(l-1)) - Hx(i,j,k-l))
                temp2(l) = temp2(l-1) + alpha(ln,l) * (Hz(i+(l-1),j,k) - Hz(i-l,j,k)) 
                enddo
                Ey(i,j,k) = Ey(i,j,k) + CEYLZ(i,j,k)*temp1(ln) + CEYLX(i,j,k)*temp2(ln)
            enddo
        enddo
    enddo

    !ソース項
    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)
            end subroutine EYFIELD



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Ez-field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EZFIELD(istep,t,Je,Ez,Hx,Hy,sigma)
    use const_para
    implicit none

    integer :: l
    integer, intent(in) :: istep
    real(8), intent(in) :: t !経過時間
!   complex(kind(0d0)),intent(in) :: sigmazz(nx,ny,nz)
    real(8),intent(in) :: Je(nstep)
    complex(kind(0d0)), intent(inout) :: Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hy(nx,ny,nz)
    real(8),intent(in)  :: sigma(nx,ny,nz)
    real(8) :: etazz(nx,ny,nz)
    real(8) :: CEZLX(nx,ny,nz)
    real(8) :: CEZLY(nx,ny,nz)
    real(8) :: alpha(ln,ln)
    real(8) :: temp1(0:ln),temp2(0:ln)
    temp1(0) = 0.0d0
    temp2(0) = 0.0d0

   !Holberg optimization scheme
    alpha(1,1) = 1.00235d0
   ! alpha(2,1:2) = (/1.14443d0,-0.04886d0/)
   ! alpha(3,1:3) = (/1.20282d0,-0.08276d0,0.00950d0/)
   ! alpha(4,1:4) = (/1.23041d0,-0.10313d0,0.02005d0,-0.00331d0/)

    
    !係数の設定
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                etazz(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmazz(i,j,k)
                CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
                CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
            enddo
        enddo
    enddo


    !波動伝播計算
    do k = 1,nz-1
        do j = 2,ny-1
            do i = 2,nx-1
                     do l=1,ln
                     temp1(l) = temp1(l-1) + alpha(ln,l) * (Hy(i+(l-1),j,k) - Hy(i-l,j,k))
                     temp2(l) = temp2(l-1) + alpha(ln,l) * (Hx(i,j+(l-1),k) - Hx(i,j-l,k))
                     enddo
                     Ez(i,j,k) = Ez(i,j,k) + CEZLX(i,j,k)*temp1(ln) + CEZLY(i,j,k)*temp2(ln)
            enddo
        enddo
    enddo

    !ソース項
!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)
            end subroutine EZFIELD
