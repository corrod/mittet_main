!!!仮想領域での電磁場の計算*****************************************************
!!!!!!!!!!!!!
!電場計算***********************************************************************
!Ex-field
subroutine EXFIELD(istep,t,Je,Ex,Hy,Hz,sigma)
    use const_para
    implicit none

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
!
    !係数の設定
    do k=1,nz
         do j=1,ny
              do i=1,nx
                    etaxx(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmaxx(i,j,k)
                    CEXLY(i,j,k) = dt * etaxx(i,j,k) / dy
                    CEXLZ(i,j,k) = - dt * etaxx(i,j,k) / dz
              enddo
         enddo
    enddo

    !電場ソースの設定
    !Je(istep) = dt * etaxx(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                  Ex(i,j,k) = Ex(i,j,k)&
                            + CEXLY(i,j,k) * (Hz(i,j,k) - Hz(i,j-1,k))&
                            + CEXLZ(i,j,k) * (Hy(i,j,k) - Hy(i,j,k-1))
            enddo
        enddo
    enddo

    !ソース項
    !Ex(x0,y0,z0) = Ex(x0,y0,z0) - Je(istep)

   !ソース位置、ソースから離れた位置でのExの磁場の時間分布
   write(9,*) 't,ex(x0+10,y0,z0)',t,real(ex(x0+10,y0,z0))
            endsubroutine EXFIELD


!Ey-field-------------------------------------------------------------
subroutine EYFIELD(istep,t,Je,Ey,Hz,Hx,sigma)
    use const_para
    implicit none

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

    !係数の設定
    do k=1,nz
        do j=1,ny
           do i=1,nx
              etayy(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmayy(i,j,k)
              CEYLZ(i,j,k) = dt * etayy(i,j,k) / dz
              CEYLX(i,j,k) = - dt * etayy(i,j,k) / dx
           enddo
        enddo
    enddo

    !電場ソースの設定
    !  Je(istep) = dt * etayy(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                Ey(i,j,k) = Ey(i,j,k)&
                            + CEYLZ(i,j,k) * (Hx(i,j,k) - Hx(i,j,k-1))&
                            + CEYLX(i,j,k) * (Hz(i,j,k) - Hz(i-1,j,k))
            enddo
        enddo
    enddo

    !ソース項
    !Ey(x0,y0,z0) = Ey(x0,y0,z0) - Je(istep)
            endsubroutine EYFIELD



!Ez-field-------------------------------------------------------------------
subroutine EZFIELD(istep,t,Je,Ez,Hx,Hy,sigma)
    use const_para
    implicit none

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

    !係数の設定
    do k=1,nz
        do j=1,ny
            do i=1,nx
                etazz(i,j,k) = 2.0d0 * omega0 * sigma(i,j,k)!sigmazz(i,j,k)
                CEZLX(i,j,k) = dt * etazz(i,j,k) / dx
                CEZLY(i,j,k) = - dt * etazz(i,j,k) / dy
            enddo
        enddo
    enddo

    !電場ソースの設定
 !   Je(istep) = dt * etazz(x0,y0,z0) * Jn(istep) /dx/dy/dz

    !波動伝播計算
    do k=1,nz-1
        do j=2,ny-1
            do i=2,nx-1
                Ez(i,j,k) = Ez(i,j,k)&
                            + CEZLX(i,j,k) * (Hy(i,j,k) - Hy(i-1,j,k))&
                            + CEZLY(i,j,k) * (Hx(i,j,k) - Hx(i,j-1,k))
            enddo
        enddo
    enddo

    !ソース項
!   Ez(x0,y0,z0) = Ez(x0,y0,z0) - Je(istep)
            endsubroutine EZFIELD
