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
    real(8), intent(in) :: myu(nx,ny,nz)
    real(8), intent(in) :: Jh(nstep)
    complex(kind(0d0)), intent(inout) :: Ex(nx,ny,nz), Ey(nx,ny,nz), Ez(nx,ny,nz)
    complex(kind(0d0)), intent(inout) :: Hx(nx,ny,nz), Hy(nx,ny,nz), Hz(nx,ny,nz)
    real(8) :: CHXLY(nx,ny,nz), CHYLZ(nx,ny,nz), CHZLX(nx,ny,nz)
    real(8) :: CHXLZ(nx,ny,nz), CHYLX(nx,ny,nz), CHZLY(nx,ny,nz)
    real(8) :: alpha(ln,ln)
    real(8) :: temp1(0:ln), temp2(0:ln)
    temp1(0) = 0.0d0
    temp2(0) = 0.0d0

    !Holberg optimization scheme
     alpha(1,1:4) = (/1.00235d0, 0.0d0, 0.0d0, 0.0d0/)
  !  alpha(2,1:4) = (/1.14443d0, -0.04886d0, 0.0d0, 0.0d0/)
  !  alpha(3,1:4) = (/1.20282d0, -0.08276d0, 0.00950d0, 0.0d0/)
  !  alpha(4,1:4) = (/1.23041d0, -0.10313d0, 0.02005d0, -0.00331d0/)


    !Hx係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nz
       do j = 1,ny
           do i = 1,nx
              CHXLY(i,j,k) = - dt / myu(i,j,k) / dy
              CHXLZ(i,j,k) = dt / myu(i,j,k) / dz
          enddo
      enddo
    enddo

    !Hx波動伝播計
    do k = 1,nz-1
        do j = 1,ny-1
            do i = 2,nx-1
                do l = 1,ln
                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ez(i,j+l,k) - Ez(i,j-(l-1),k))
                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ey(i,j,k+l) - Ey(i,j,k-(l-1)))
                enddo
                Hx(i,j,k) = Hx(i,j,k) + CHXLY(i,j,k)*temp1(ln) + CHXLZ(i,j,k)*temp2(ln)
            enddo
        enddo
    enddo

    !ソース項
!   Hx(x0,y0,z0) = Hx(x0,y0,z0) - Jh(istep)


	!Hy係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
             CHYLZ(i,j,k) = - dt / myu(i,j,k) / dz
             CHYLX(i,j,k) = dt / myu(i,j,k) / dx
          enddo
      enddo
    enddo

    !Hy波動伝播計算
    do k = 1,nz-1
        do j = 2,ny-1
            do i = 1,nx-1
                do l = 1,ln
                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ex(i,j,k+l) - Ex(i,j,k-(l-1)))
                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ez(i+l,j,k) - Ez(i-(l-1),j,k))
                enddo
                Hy(i,j,k) = Hy(i,j,k) + CHYLZ(i,j,k)*temp1(ln) + CHYLX(i,j,k)*temp2(ln)
            enddo
        enddo
    enddo

    !ソース項
	!   Hy(x0,y0,z0) = Hy(x0,y0,z0) - Jh(istep)


	!Hz係数の設定!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do k = 1,nz
       do j = 1,ny
         do i = 1,nx
            CHZLX(i,j,k) = - dt / myu(i,j,k) / dz
            CHZLY(i,j,k) = dt / myu(i,j,k) / dx
        enddo
      enddo
    enddo

    !Hz波動伝播計算
    do k = 2,nz-1
        do j = 1,ny-1
            do i = 1,nx-1
                do l = 1,ln
                    temp1(l) = temp1(l-1) + alpha(ln,l) * (Ey(i+l,j,k) - Ey(i-(l-1),j,k))
                    temp2(l) = temp2(l-1) + alpha(ln,l) * (Ex(i,j+l,k) - Ex(i,j-(l-1),k))
                enddo
                Hz(i,j,k) = Hz(i,j,k) + CHZLX(i,j,k)*temp1(ln) + CHZLY(i,j,k)*temp2(ln) 
             enddo           
        enddo
    enddo

    !ソース項
    Hz(x0,y0,z0) = Hz(x0,y0,z0) - Jh(istep)
            end subroutine Hfield